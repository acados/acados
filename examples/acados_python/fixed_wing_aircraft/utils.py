import csv
import json

import casadi as ca
import numpy as np


def load_trajectory(filename: str, dt: float, number_of_orbits: float = 1, drop_t_states: bool = True) -> tuple[np.ndarray, np.ndarray, float]:
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)
        data = np.array([list(map(float, row)) for row in reader])

    time_col = data[:, headers.index('X_t_0')]

    x_idx  = [i for i, h in enumerate(headers) if h.startswith('X_')]
    u_idx  = [i for i, h in enumerate(headers) if h.startswith('U_')]

    if drop_t_states:
        x_idx  = [i for i in  x_idx if not any(p in headers[i] for p in ('_t_', '_dt_'))]

    orbit_duration = time_col[-1]
    total_duration = number_of_orbits * orbit_duration

    n_steps = int(np.ceil(total_duration / dt))  
    t_absolute = np.arange(n_steps + 1) * dt

    # the trajectory is periodic, so fold the absolute times back into [0, T)
    t_query = t_absolute % orbit_duration

    def resample(col_idx: list[int]) -> np.ndarray:
        columns = [np.interp(t_query, time_col, data[:, i]) for i in col_idx]
        return np.column_stack(columns)

    return resample(x_idx), resample(u_idx), orbit_duration


def wind_profile(phases: list[np.ndarray], gust_duration: float, t: np.ndarray) -> np.ndarray:
    """Build a piecewise-constant wind with a 1-cos gust at every switch.

    The phases hold one wind vector each and are spread evenly over `t`. The
    wind starts at rest and gusts into the first phase, then into each following
    one at its interval boundary, blending over `gust_duration` as
    0.5 * (1 - cos), which starts and ends with zero slope.

    Parameters
    ----------
    phases : list of ndarray (3,)
        Wind vectors in NED, one per interval.
    gust_duration : float
        Length of the transition at each switch, in the units of `t`. Zero makes
        the wind start at the first phase and switch instantaneously.
    t : ndarray (N,)
        Time grid to evaluate the profile on.

    Returns
    -------
    wind : ndarray (N, 3)
    """
    phases = np.asarray(phases, dtype=float)
    t = np.asarray(t, dtype=float)

    # interval boundaries, evenly spaced over the grid
    edges = np.linspace(t[0], t[-1], len(phases) + 1)

    # piecewise-constant base: the interval each sample falls into
    interval = np.clip(np.searchsorted(edges, t, side='right') - 1, 0, len(phases) - 1)
    wind = phases[interval]

    # each gust starts at its boundary and leaves the phase before it, the first
    # one starting from rest; without a gust the phases switch instantaneously
    previous = np.vstack([np.zeros_like(phases[0]), phases[:-1]])
    for k in range(len(phases) if gust_duration > 0 else 0):
        s = (t - edges[k]) / gust_duration
        gust = (s >= 0.0) & (s < 1.0)
        blend = 0.5 * (1.0 - np.cos(np.pi * s[gust]))
        wind[gust] = previous[k] + (phases[k] - previous[k]) * blend[:, np.newaxis]

    return wind


def buildForcesanMomentsFunction(coeff_file: str) -> ca.Function:
    """

    Parameters
    ----------
    coeff_file: path to the coefficients file RELATIVE TO THE BASE DIRECTORY

    Returns
    -------

    """

    with open(coeff_file, 'r') as f:
        data = json.load(f)

    # create symbolic variables for the inputs: u,v,w,p,q,r,Thrust_L,Thrust_R,Elevator,Aileron_L,Aileron_R
    velocities = ca.SX.sym('velocities', 3)  # [u,v,w]
    rates = ca.SX.sym('rates', 3)  # [p,q,r]
    thrust = ca.SX.sym('thrust', 2)  # [Thrust_L, Thrust_R]
    control_surfaces = ca.SX.sym('control_surfaces', 2)  # [Elevator, Aileron_LR]

    # 1. Define Basic Symbolic Inputs (matching your FEATURES list)
    U,V,W = velocities[0], velocities[1], velocities[2]
    p,q,r = rates[0], rates[1], rates[2]
    TL, TR = thrust[0], thrust[1]
    Elev = control_surfaces[0]
    AilL = control_surfaces[1]

    # Store base features in a dictionary for easy interaction generation
    base_map = {
        'U': U, 'V': V, 'W': W,
        'p': p, 'q': q, 'r': r,
        'Thrust_L': TL, 'Thrust_R': TR,
        'Elevator': Elev, 'Aileron_L': AilL
    }

    # 2. Replicate the "Feature Engineering" logic
    # This dictionary maps the JSON string keys to CasADi expressions
    feat_expr = {**base_map}

    # Add Squared terms: f^2
    for name, expr in base_map.items():
        feat_expr[f"{name}^2"] = expr ** 2

    # Add Interaction terms: f*g (using same sorting logic as Polars snippet)
    names = sorted(base_map.keys())
    for i, f in enumerate(names):
        for g in names[i + 1:]:
            feat_expr[f"{f}{g}"] = base_map[f] * base_map[g]

    # Add the special Aerodynamic features
    alpha = W / U
    feat_expr["alpha"] = alpha
    feat_expr["alpha^2"] = alpha ** 2
    feat_expr["alpha^3"] = alpha ** 3
    feat_expr["W3/U"] = (W ** 3) / U
    feat_expr["U2_Elevator"] = (U ** 2) * Elev
    feat_expr["U2_Aileron_L"] = (U ** 2) * AilL
    feat_expr["offset"] = 1.0

    # 3. Build the Forces and Moments
    results = {}
    for target in ['Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']:
        coeffs = data['coefficients'][target]
        expr = 0
        for name, val in coeffs.items():
            if val != 0:  # Sparse optimization: skip zero coefficients
                expr += val * feat_expr[name]
        results[target] = expr

    # 4. Create the CasADi Function
    # Input order: [U, V, W, p, q, r, TL, TR, Elev, AilR]
    f_aero = ca.Function('aero_model',
                         [velocities, rates, thrust, control_surfaces],
                         [ca.vertcat(results['Fx'], results['Fy'], results['Fz']),
                          ca.vertcat(results['Mx'], results['My'], results['Mz'])],
                         ['velocities', 'rates', 'thrust', 'control_surfaces'],
                         ['forces', 'moments'])
    return f_aero


def rotm2eul(R):
    """
    Extract yaw, pitch, and roll angles from rotation matrix/matrices.

    Computes Euler angles from one or more 3x3 rotation matrices using
    the ZYX convention (yaw-pitch-roll, intrinsic rotations).

    Parameters
    ----------
    R : ndarray
        Rotation matrix of shape (3, 3) for a single rotation, or
        (N, 3, 3) for N rotation matrices.

    Returns
    -------
    yaw : float or ndarray
        Rotation angle around the z-axis in radians. Scalar if R is (3, 3),
        array of shape (N,) if R is (N, 3, 3).
    pitch : float or ndarray
        Rotation angle around the y-axis in radians. Scalar if R is (3, 3),
        array of shape (N,) if R is (N, 3, 3).
    roll : float or ndarray
        Rotation angle around the x-axis in radians. Scalar if R is (3, 3),
        array of shape (N,) if R is (N, 3, 3).
    """
    match R.shape:
        case (3, 3):
            roll  = np.arctan2(R[2, 1], R[2, 2])
            pitch = np.arcsin(np.clip(-R[2, 0], -1.0, 1.0))
            yaw   = np.arctan2(R[1, 0], R[0, 0])
        case (_, 3, 3):
            roll  = np.arctan2(R[:, 2, 1], R[:, 2, 2])
            pitch = np.arcsin(np.clip(-R[:, 2, 0], -1.0, 1.0))
            yaw   = np.arctan2(R[:, 1, 0], R[:, 0, 0])
        case _:
            raise ValueError(f"Expected R to have shape (3, 3) or (N, 3, 3), got {R.shape}")
    return yaw, pitch, roll


def inertial_to_body(v_inertial, R):
    """
    Transform velocity from inertial to body frame.

    Accepts a single vector/matrix pair of shapes (3,) and (3, 3), or N samples
    of shapes (N, 3) and (N, 3, 3). The result has the shape of `v_inertial`.
    """
    match (v_inertial.shape, R.shape):
        case ((3,), (3, 3)):
            return R.T @ v_inertial
        case ((n, 3), (m, 3, 3)) if n == m:
            return np.einsum('kji,kj->ki', R, v_inertial)
        case ((_, 3), (_, 3, 3)):
            raise ValueError(
                f"Dimension mismatch: v_inertial has {v_inertial.shape[0]} samples, "
                f"R has {R.shape[0]} samples"
            )
        case _:
            raise ValueError(
                f"Expected shapes (3,) and (3, 3), or (N, 3) and (N, 3, 3), "
                f"got {v_inertial.shape} and {R.shape}"
            )

def angle_of_attack(v_inertial, R):
    """
    Compute the angle of attack from inertial velocity and rotation matrix.

    Transforms the velocity vector from inertial frame to body frame and
    computes the angle of attack as the angle between the velocity vector
    and the body x-axis in the xz-plane.

    Parameters
    ----------
    v_inertial : ndarray
        Velocity vector in inertial frame. Shape (3,) for a single vector,
        or (N, 3) for N velocity vectors.
    R : ndarray
        Rotation matrix from body to inertial frame. Shape (3, 3) for a
        single rotation, or (N, 3, 3) for N rotation matrices.

    Returns
    -------
    alpha : float or ndarray
        Angle of attack in radians. Scalar if inputs are single vector/matrix,
        array of shape (N,) for batched inputs.
    """
    v_body = inertial_to_body(v_inertial, R)
    return np.atan2(v_body[..., 2], v_body[..., 0])


def sideslip_angle(v_inertial, R):
    """
    Compute the sidesip angle from inertial velocity and rotation matrix.

    Transforms the velocity vector from inertial frame to body frame and
    computes the sideslip angle as the angle between the velocity vector
    and the body x-axis in the xy-plane.

    Parameters
    ----------
    v_inertial : ndarray
        Velocity vector in inertial frame. Shape (3,) for a single vector,
        or (N, 3) for N velocity vectors.
    R : ndarray
        Rotation matrix from body to inertial frame. Shape (3, 3) for a
        single rotation, or (N, 3, 3) for N rotation matrices.

    Returns
    -------
    alpha : float or ndarray
        Sideslip angle in radians. Scalar if inputs are single vector/matrix,
        array of shape (N,) for batched inputs.
    """
    v_body = inertial_to_body(v_inertial, R)
    # approx. for small aoa; exact is atan2(v, hypot(u, w)), matches the model
    return np.atan2(v_body[..., 1], v_body[..., 0])
