#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import blended_transform_factory
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from utils import angle_of_attack, inertial_to_body, rotm2eul, sideslip_angle

from acados_template import latexify_plot


latexify_plot()

# Transformation matrix from NED to ENU for plotting
R_NED_TO_ENU = np.array([[0, 1, 0],
                         [1, 0, 0],
                         [0, 0, -1]])

# Constraint bounds
OMEGA_BOUNDS = (-1.0, 1.0)                      # [rad/s]
THRUST_BOUNDS = (0.0, 20.0)                     # [N]
ELEVATOR_BOUNDS = (-32.0, 32.0)                 # [deg]
AILERON_BOUNDS = (-25.0, 25.0)                  # [deg]
DTHRUST_BOUNDS = (-5.0, 5.0)                    # [N/s]
DFLAPS_BOUNDS = (-19.0, 19.0)                   # [deg/s]
AOA_BOUNDS = (-6.0, 12.0)                       # [deg]
SIDESLIP_BOUNDS = (-10.0, 10.0)                 # [deg]

# Initial 3D view angles [deg]
VIEW_ELEVATION = 30
VIEW_AZIMUTH = -65

# line styles
TRUE_STYLE = dict(linestyle='-')                                # simulated truth
REFERENCE_STYLE = dict(linestyle='--', alpha=0.5)               # setpoint
ESTIMATE_STYLE = dict(linestyle=':', linewidth=1.5)             # EKF estimate
PREDICTION_STYLE = dict(linestyle='-', linewidth=1, alpha=0.7)  # MPC horizon
BAND_STYLE = dict(alpha=0.2, linewidth=0)                       # +-n_sigma
BOUND_STYLE = dict(color='k', linestyle='--', alpha=0.5)        # constraint


# ==================== Data helpers ====================

def _rotation_matrices(R_flat):
    """Read the flattened R block of a state trajectory as (N, 3, 3) matrices.

    The state holds R column-major, so the samples are transposed after reshaping.
    """
    return np.asarray(R_flat).reshape(-1, 3, 3).transpose(0, 2, 1)


def _unpack_states(X):
    """Split a (N, 22) state trajectory into its named components.

    The rotation matrices are returned both as (N, 3, 3) and as the euler
    angles (in radians) they correspond to.
    """
    X = np.asarray(X)
    R = _rotation_matrices(X[:, 6:15])
    yaw, pitch, roll = rotm2eul(R)

    return dict(
        r=X[:, 0:3],
        dr=X[:, 3:6],
        R=R,
        omega=X[:, 15:18],
        thrust=X[:, 18:20],
        flaps=X[:, 20:22],
        roll=roll, pitch=pitch, yaw=yaw,
    )


def _to_enu(v):
    """Convert a (N, 3) array of NED vectors to ENU."""
    return (R_NED_TO_ENU @ np.asarray(v).T).T


def _add_bounds(ax, bounds):
    """Draw horizontal constraint lines."""
    for bound in bounds:
        ax.axhline(bound, **BOUND_STYLE)


def _ylim(*arrays, margin_frac=0.1):
    """Common y-limits over all given arrays, with a relative margin.

    Missing and empty arrays are skipped, so a panel without a prediction,
    reference or constraint bounds can pass them anyway.
    """
    values = np.concatenate([np.asarray(a).ravel() for a in arrays
                             if a is not None and np.size(a) > 0])

    lower, upper = values.min(), values.max()
    span = max(upper - lower, 0.1)

    return lower - margin_frac * span, upper + margin_frac * span


def _body_velocity(X, wind):
    """Velocity relative to the air in the body frame, (N, 3)."""
    R = _rotation_matrices(X[:, 6:15])
    return inertial_to_body(X[:, 3:6] - wind, R)


def _draw_panel(ax, t, labels, ylabel, series, bounds=(), scale=1.0, step=False):
    """Draw one panel and return one line group per entry in `series`.

    `series` holds (data, style) pairs in draw order. A `data` of None creates
    empty lines for the caller to fill in per animation frame. The labels go on
    the first entry, so the legend swatches belong to the panel's primary line.
    """
    draw = ax.step if step else ax.plot
    step_kw = dict(where='post') if step else {}

    groups = []
    for position, (data, style) in enumerate(series):
        group = []
        for k in range(len(labels)):
            label = labels[k] if position == 0 else None
            if data is None:
                line, = draw([], [], color=f'C{k}', label=label, **style, **step_kw)
            else:
                line, = draw(t, data[:, k] * scale, color=f'C{k}', label=label,
                             **style, **step_kw)
            group.append(line)
        groups.append(group)

    _add_bounds(ax, bounds)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=len(labels), loc='upper right')

    return groups


def _plot_series(ax, t, labels, ylabel, true=None, reference=None, estimate=None,
                 std=None, bounds=(), scale=1.0, step=False, n_sigma=1.0):
    """Draw one static panel, each series in the style of the role it plays.

    One colour per label; `std` shades a +-`n_sigma` band around the estimate.
    """
    # the band first, so it stays behind every line
    if std is not None:
        for k in range(len(labels)):
            band = n_sigma * std[:, k] * scale
            ax.fill_between(t, estimate[:, k] * scale - band, estimate[:, k] * scale + band,
                            color=f'C{k}', **BAND_STYLE)

    roles = ((true, TRUE_STYLE), (estimate, ESTIMATE_STYLE), (reference, REFERENCE_STYLE))
    series = [(data, style) for data, style in roles if data is not None]

    _draw_panel(ax, t, labels, ylabel, series, bounds=bounds, scale=scale, step=step)


def _style_legend(fig, roles, **kwargs):
    """Name the line styles once for the whole figure, in neutral grey.

    `roles` are (label, style) pairs; a style without a linestyle is taken as a
    shaded band and gets a patch handle instead of a line. Grey stands in for
    the component colours, unless the style prescribes one of its own.
    """
    handles = [Patch(**{'color': 'gray', **style}) if 'linestyle' not in style
               else Line2D([], [], **{'color': 'gray', **style})
               for _, style in roles]

    fig.legend(handles, [label for label, _ in roles], loc='lower center',
               ncol=len(roles), frameon=False, **kwargs)


def _series_roles(reference=False, estimate=False, bounds=False, n_sigma=1.0):
    """The `_style_legend` roles for a plot of true values and its companions."""
    roles = [('true', TRUE_STYLE)]
    if reference:
        roles.append(('reference', REFERENCE_STYLE))
    if estimate:
        roles += [('estimate', ESTIMATE_STYLE), (f'$\\pm{n_sigma:g}\\sigma$', BAND_STYLE)]
    if bounds:
        roles.append(('bound', BOUND_STYLE))

    return roles


# ==================== 3D drawing helpers ====================

# Airplane silhouette in body coordinates, one column per vertex.
BODY_VERTICES = np.array([
    [0.35, 0, 0],
    [0.4, 0, -0.1],
    [0, 0, -0.15],
    [-0.5, 0, -0.1],
    [-0.6, 0, -0.2],
    [-0.7, 0, -0.2],
    [-0.6, 0, 0]
]).T + np.array([[0, 0, 0.05]]).T

WINGS_VERTICES = np.array([
    [0.25, 0, 0],
    [0.2, 0.5, 0],
    [0, 0.5, 0],
    [0, -0.5, 0],
    [0.2, -0.5, 0]
]).T


def _draw_plane(ax, position_NED, R_NED, scale=1.0, facecolor='lightblue',
                alpha=0.8, draw_axes=True):
    """Draw a 3D airplane, returning the artists that were added."""
    position_NED = np.asarray(position_NED).flatten().reshape((3, 1))
    R_NED = np.asarray(R_NED).reshape(3, 3)

    position_ENU = R_NED_TO_ENU @ position_NED

    # Transform to world frame, then to ENU for plotting
    body_vertices_ENU = R_NED_TO_ENU @ (R_NED @ BODY_VERTICES * scale + position_NED)
    wings_vertices_ENU = R_NED_TO_ENU @ (R_NED @ WINGS_VERTICES * scale + position_NED)

    artists = []
    for vertices_ENU in (body_vertices_ENU, wings_vertices_ENU):
        collection = Poly3DCollection([vertices_ENU.T], edgecolor='black',
                                      facecolor=facecolor, alpha=alpha)
        ax.add_collection3d(collection)
        artists.append(collection)

    # Draw body coordinate axes
    if draw_axes:
        R_plot = R_NED_TO_ENU @ R_NED
        for i, c in enumerate(['C0', 'C1', 'C2']):
            quiver = ax.quiver(*position_ENU.flatten(), *(scale * R_plot[:, i]),
                               color=c, arrow_length_ratio=0.2)
            artists.append(quiver)

    return artists


def _trajectory_bounds(positions_ENU):
    """Cubic bounds around the trajectory: (centre, edge length).

    The three axes share one range, so distances stay comparable in every
    direction.
    """
    positions_ENU = np.atleast_2d(positions_ENU)
    lower = positions_ENU.min(axis=0)
    upper = positions_ENU.max(axis=0)

    return (lower + upper) / 2, max((upper - lower).max(), 1.0) * 1.1


def _setup_3d_axes(ax, center, max_range):
    """Configure 3D axes with cubic, trajectory-centred bounds."""
    ax.set_xlim(center[0] - max_range / 2, center[0] + max_range / 2)
    ax.set_ylim(center[1] - max_range / 2, center[1] + max_range / 2)
    ax.set_zlim(center[2] - max_range / 2, center[2] + max_range / 2)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel('Y [m] (East)')
    ax.set_ylabel('X [m] (North)')
    ax.set_zlabel('-Z [m] (Up)')
    ax.grid(True, alpha=0.25)
    ax.view_init(elev=VIEW_ELEVATION, azim=VIEW_AZIMUTH)


def _draw_wind_arrow(ax, wind_NED, center, max_range, arrow_scale=0.15, label='wind'):
    """Draw a wind arrow indicator in the lower corner of the viewport.

    `label` names what the arrow shows, since the static plot draws the mean over
    the run while the animation draws the wind of the current frame.
    """
    artists = []

    v_w_ENU = R_NED_TO_ENU @ np.asarray(wind_NED).flatten()
    wind_magnitude = np.linalg.norm(v_w_ENU)

    if wind_magnitude < 0.1:
        return artists

    arrow_origin = center + max_range * np.array([-0.45, -0.45, -0.45])
    arrow_vector = v_w_ENU * arrow_scale * max_range / max(wind_magnitude, 1.0)

    artists.append(ax.quiver(
        *arrow_origin, *arrow_vector,
        color='darkblue', arrow_length_ratio=0.3, linewidth=2.5, alpha=0.8
    ))

    label_pos = arrow_origin + np.array([0, 0, max_range * 0.07])
    artists.append(ax.text(
        *label_pos, f'{label}: {wind_magnitude:.1f} m/s',
        color='darkblue', fontsize=10, ha='center', va='bottom'
    ))

    return artists


# ==================== Plots ====================

def plot_states(X, t, x_ref=None, wind=0.0, show_bounds=True, figsize=(7, 12), show=False):
    """Plot the state trajectories and the aerodynamic angles.

    Thrust and flaps are states too, but live in `plot_controls` next to their
    rates. Angle of attack and sideslip are not states, but fit here well - they are what
    the nonlinear constraints act on.

    X (N, 22) states, t (>=N,) time, x_ref (>=N, 22) reference drawn dashed,
    wind (>=N, 3) in NED - the aerodynamic angles are relative to the air, and
    the reference is a no-wind trajectory so its angles use still air.
    """
    states = _unpack_states(X)
    n_samples = len(X)
    t = np.asarray(t)[:n_samples]
    wind = np.asarray(wind)[:n_samples] if np.ndim(wind) else wind

    ref = _unpack_states(np.asarray(x_ref)[:n_samples]) if x_ref is not None else None

    def aero_angles(s, w):
        """Angle of attack and sideslip as (N, 1) columns, relative to the air."""
        v = s['dr'] - w
        return angle_of_attack(v, s['R'])[:, None], sideslip_angle(v, s['R'])[:, None]

    aoa, sideslip = aero_angles(states, wind)
    ref_aoa, ref_sideslip = aero_angles(ref, 0.0) if ref is not None else (None, None)

    euler = lambda s: np.column_stack([s['roll'], s['pitch']])

    fig, axes = plt.subplots(nrows=6, ncols=1, num='States', figsize=figsize, sharex=True)

    def plot_row(ax, true, reference, labels, ylabel, bounds=(), scale=1.0):
        _plot_series(ax, t, labels, ylabel, true=true, reference=reference,
                     bounds=bounds if show_bounds else (), scale=scale)

    def both(key):
        """A state block and its reference, or None where there is no reference."""
        return states[key], None if ref is None else ref[key]

    plot_row(axes[0], *both('r'), ['x', 'y', 'z'], 'position [m]')
    plot_row(axes[1], *both('dr'), ['vx', 'vy', 'vz'], 'velocity [m/s]')
    plot_row(axes[2], *both('omega'), ['p', 'q', 'r'], 'angular velocity [rad/s]',
             bounds=OMEGA_BOUNDS)
    plot_row(axes[3], euler(states), None if ref is None else euler(ref),
             ['roll', 'pitch'], 'euler angles [deg]', scale=np.rad2deg(1))
    plot_row(axes[4], aoa, ref_aoa, [r'$\alpha$'], 'angle of attack [deg]',
             bounds=AOA_BOUNDS, scale=np.rad2deg(1))
    plot_row(axes[5], sideslip, ref_sideslip, [r'$\beta$'], 'sideslip angle [deg]',
             bounds=SIDESLIP_BOUNDS, scale=np.rad2deg(1))

    axes[-1].set_xlabel('time [s]')

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    _style_legend(fig, _series_roles(reference=x_ref is not None, bounds=show_bounds))
    if show:
        plt.show()

    return fig, axes


def plot_estimations(X, X_est, t, P, wind, wind_est, n_sigma=1.0,
                     figsize=(16, 9), show=False):
    """Plot the EKF estimate with its standard deviation against the truth.

    Truth solid, estimate dotted with its +-`n_sigma` band. Every series shown
    is an estimator state, so its band comes off the diagonal of `P`.

    X, X_est (N, 22) true and estimated states, t (>=N,) time,
    P (N, 25, 25) covariance of [states, wind], wind (>=N, 3) and
    wind_est (N, 3) true and estimated wind in NED.
    """
    n_samples = min(len(X), len(X_est))
    t = np.asarray(t)[:n_samples]

    # `P` is over the augmented estimator state [states (22), wind (3)]
    true = np.column_stack([np.asarray(X)[:n_samples], np.asarray(wind)[:n_samples]])
    est = np.column_stack([np.asarray(X_est)[:n_samples], np.asarray(wind_est)[:n_samples]])
    std = np.sqrt(np.maximum(np.diagonal(np.asarray(P)[:n_samples], axis1=1, axis2=2), 0.0))

    fig, axes = plt.subplots(nrows=3, ncols=2, num='Estimation', figsize=figsize, sharex=True)

    def plot_panel(ax, columns, labels, ylabel, scale=1.0):
        _plot_series(ax, t, labels, ylabel, true=true[:, columns],
                     estimate=est[:, columns], std=std[:, columns],
                     scale=scale, n_sigma=n_sigma)

    plot_panel(axes[0, 0], [0, 1, 2], ['x', 'y', 'z'], 'position [m]')
    plot_panel(axes[0, 1], [18, 19], ['left', 'right'], 'thrust [N]')

    plot_panel(axes[1, 0], [3, 4, 5], ['vx', 'vy', 'vz'], 'velocity [m/s]')
    plot_panel(axes[1, 1], [20, 21], ['elevator', 'aileron'], 'flap deflection [deg]',
               scale=np.rad2deg(1))

    plot_panel(axes[2, 0], [15, 16, 17], ['p', 'q', 'r'], 'angular velocity [rad/s]')
    plot_panel(axes[2, 1], [22, 23, 24], ['$w_x$', '$w_y$', '$w_z$'], 'wind [m/s]')

    for ax in axes[-1]:
        ax.set_xlabel('time [s]')
    fig.suptitle('EKF estimate vs. truth')

    fig.tight_layout(rect=(0, 0.04, 1, 1))
    _style_legend(fig, _series_roles(estimate=True, n_sigma=n_sigma))
    if show:
        plt.show()

    return fig, axes


def plot_controls(U, X, t, u_ref=None, x_ref=None, show_bounds=True,
                  figsize=(7, 10), show=False):
    """Plot each actuator next to the rate that drives it.

    Thrust and flap deflection are model states, the controls are their rates,
    so each pair shares a row: thrust over dthrust, flaps over dflaps.

    U (M, 4) controls [dthrust_left, dthrust_right, delevator, daileron],
    X (N, 22) states, t (>=N,) time, u_ref (>=M, 4) and x_ref (>=N, 22)
    references drawn dashed. State rows span N samples, rate rows M.
    """
    U = np.asarray(U)
    n_u, n_x = len(U), len(X)
    t = np.asarray(t)
    t_u, t_x = t[:n_u], t[:n_x]

    states = _unpack_states(X)
    ref = _unpack_states(np.asarray(x_ref)[:n_x]) if x_ref is not None else None
    u_ref = np.asarray(u_ref)[:n_u] if u_ref is not None else None

    fig, axes = plt.subplots(nrows=4, ncols=1, num='Controls', figsize=figsize, sharex=True)

    def plot_state_row(ax, key, labels, ylabel, bounds, scale=1.0):
        _plot_series(ax, t_x, labels, ylabel, true=states[key],
                     reference=None if ref is None else ref[key],
                     bounds=bounds if show_bounds else (), scale=scale)

    def plot_rate_row(ax, columns, labels, ylabel, bounds, scale=1.0):
        # the MPC applies piecewise constant inputs
        _plot_series(ax, t_u, labels, ylabel, true=U[:, columns],
                     reference=None if u_ref is None else u_ref[:, columns],
                     bounds=bounds if show_bounds else (), scale=scale, step=True)

    plot_state_row(axes[0], 'thrust', ['left', 'right'], 'thrust [N]', THRUST_BOUNDS)
    plot_rate_row(axes[1], [0, 1], ['left', 'right'], 'dthrust [N/s]', DTHRUST_BOUNDS)
    plot_state_row(axes[2], 'flaps', ['elevator', 'left=-right aileron'],
                   'flap deflection [deg]', (*ELEVATOR_BOUNDS, *AILERON_BOUNDS),
                   scale=np.rad2deg(1))
    plot_rate_row(axes[3], [2, 3], ['elevator', 'left=-right aileron'], 'dflaps [deg/s]',
                  DFLAPS_BOUNDS, scale=np.rad2deg(1))

    axes[-1].set_xlabel('time [s]')

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    _style_legend(fig, _series_roles(reference=u_ref is not None or ref is not None,
                                     bounds=show_bounds))
    if show:
        plt.show()

    return fig, axes


def plot_3d_static(X, x_ref=None, wind=None, num_planes=9, scale=None,
                   figsize=(12, 10), show=False):
    """Plot the flown trajectory in 3D, with the plane drawn along it.

    X (N, 22) states, x_ref (>=N, 22) reference drawn dashed, wind (>=N, 3) in
    NED, shown as a corner arrow of its mean (a near zero wind is skipped).
    `num_planes` snapshots are drawn, `scale` sizing them (default 10% of the
    viewport).
    """
    states = _unpack_states(X)
    positions_ENU = _to_enu(states['r'])
    n_samples = len(positions_ENU)

    fig = plt.figure('3D Trajectory (Static)', figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    center, max_range = _trajectory_bounds(positions_ENU)
    _setup_3d_axes(ax, center, max_range)

    if scale is None:
        scale = 0.1 * max_range

    # Trajectory and reference
    ax.plot(*positions_ENU.T, color='C0', linewidth=2, label='trajectory')
    if x_ref is not None:
        reference_ENU = _to_enu(np.asarray(x_ref)[:n_samples, 0:3])
        ax.plot(*reference_ENU.T, color='C0', linestyle='--', alpha=0.4,
                linewidth=1.5, label='reference')

    if wind is not None:
        _draw_wind_arrow(ax, np.mean(np.asarray(wind)[:n_samples], axis=0),
                         center, max_range, label='wind avg.')

    # Airplane snapshots, fading in over time. endpoint=False leaves the final
    # state undrawn, so on a closed orbit it does not cover the red starting one
    plane_indices = np.linspace(0, n_samples - 1, num_planes, dtype=int, endpoint=False)
    for i, idx in enumerate(plane_indices):
        _draw_plane(
            ax, states['r'][idx], states['R'][idx], scale=scale,
            facecolor='red' if i == 0 else 'lightblue',
            alpha=0.3 + 0.7 * (i / len(plane_indices)),
            draw_axes=(i == 0),  # only on the first snapshot
        )

    ax.legend(loc='upper left')
    fig.tight_layout()
    if show:
        plt.show()

    return fig, ax


def plot_3d_animated(X, U, t, X_pred, U_pred, x_ref, u_ref, wind, wind_est, P,
                     n_sigma=1.0, scale=None, figsize=(16, 9), dpi=120,
                     time_scaling=1.0, frame_step=1, output_path=None):
    """
    Animate the closed loop: the 3D view on the left half, four stacked panels
    on the right (position, body velocity, wind estimate, dthrust).

    Each panel splits at the current time: two thirds past, one third the MPC
    prediction of that frame. The window is locked to the horizon, so it fills
    up and then scrolls along, keeping the 2:1 split. The prediction starts from
    the estimator state, so the offset where it meets the past is the estimation
    error.

    X (N, 22) states and U (M, 4) controls - the animation runs over min(N, M)
    frames, which also skips the trailing prediction slot `main.py` never
    writes; t (>=N,) time, its first interval setting the sampling time;
    X_pred (M, 22, N_horizon+1) and U_pred (M, 4, N_horizon) the predictions per
    step; x_ref (>=N+N_horizon, 22) and u_ref reference, drawn dashed across past
    and future alike; wind (>=N, 3) true wind, wind_est (N, 3) its estimate and
    P (N, 25, 25) the covariance behind the wind band.

    Returns (anim, fig), or (None, None) if written to `output_path`. The
    animation must stay referenced, or it is garbage collected before it plays.
    """
    X = np.asarray(X)
    U = np.asarray(U)
    t = np.asarray(t)

    n_frames = min(len(X), len(U))
    dt = float(t[1] - t[0])
    N_horizon = np.asarray(U_pred).shape[2]
    horizon = N_horizon * dt

    t_sim = t[:n_frames]

    # ===== Flown trajectory =====
    states = _unpack_states(X[:n_frames])
    positions_ENU = _to_enu(states['r'])
    wind_NED = np.asarray(wind, dtype=float)[:n_frames]
    v_body = _body_velocity(X[:n_frames], wind_NED)

    # ===== Reference, known over the past and the horizon alike =====
    ref = np.asarray(x_ref)[:n_frames + N_horizon + 1]
    ref_controls = np.asarray(u_ref)[:n_frames + N_horizon + 1]
    ref_positions_ENU = _to_enu(ref[:, 0:3])
    # the reference is a nominal no-wind trajectory, so its body velocity is
    # taken relative to still air
    ref_v_body = _body_velocity(ref, 0.0)

    # ===== Wind estimate and its standard deviation =====
    wind_estimate = np.asarray(wind_est, dtype=float)[:n_frames]
    variance = np.diagonal(np.asarray(P)[:n_frames], axis1=1, axis2=2)[:, -3:]
    wind_std = np.sqrt(np.maximum(variance, 0.0))

    # ===== Predictions =====
    pred_states = np.asarray(X_pred)[:n_frames].transpose(0, 2, 1)   # (n, N_h+1, 22)
    pred_positions = pred_states[:, :, 0:3]
    pred_dthrust = np.asarray(U_pred)[:n_frames, 0:2, :]             # (n, 2, N_h)

    # the MPC assumes the current wind estimate over the whole horizon
    held_wind = np.repeat(wind_estimate[:, np.newaxis, :], N_horizon + 1, axis=1)
    pred_v_body = _body_velocity(pred_states.reshape(-1, X.shape[1]),
                                 held_wind.reshape(-1, 3)).reshape(n_frames, -1, 3)

    # ===== Figure: 3D on the left half, panels on the right =====
    fig = plt.figure('3D Animation', figsize=figsize, dpi=dpi)
    grid = GridSpec(1, 2, figure=fig, width_ratios=[1, 1],
                    left=0.02, right=0.95, top=0.93, bottom=0.11, wspace=0.12)

    ax3d = fig.add_subplot(grid[0, 0], projection='3d')
    panel_grid = grid[0, 1].subgridspec(4, 1, hspace=0.15)
    panels = [fig.add_subplot(panel_grid[k, 0]) for k in range(4)]
    ax_position, ax_velocity, ax_wind, ax_dthrust = panels

    for ax in panels[1:]:
        ax.sharex(panels[0])
    for ax in panels[:-1]:
        ax.tick_params(labelbottom=False)
    ax_dthrust.set_xlabel('time [s]')

    # while the past fills up the window reaches back before the start of the
    # run, so those ticks are left unlabelled
    start = t_sim[0]
    ax_dthrust.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: '' if value < start else f'{value:g}'))

    # ===== Static 3D content =====
    center, max_range = _trajectory_bounds(np.vstack([positions_ENU, ref_positions_ENU]))
    _setup_3d_axes(ax3d, center, max_range)

    if scale is None:
        scale = 0.1 * max_range

    ax3d.plot(*ref_positions_ENU.T, color='C0', linewidth=1.5,
              label='reference', **REFERENCE_STYLE)
    trail, = ax3d.plot([], [], [], color='C0', linewidth=2,
                       label='trajectory', **TRUE_STYLE)
    ax3d.legend(loc='upper left')

    plane_artists = _draw_plane(ax3d, states['r'][0], states['R'][0], scale=scale)
    wind_artists = []

    # ===== Static panel content =====
    def setup_panel(ax, labels, ylabel, ylim_data, pred_data=None, ref_data=None,
                    bounds=(), step=False, roles=(TRUE_STYLE,)):
        """Draw the reference, bounds and decoration, return the animated lines.

        One animated line group per style in `roles` (the panel's own series,
        which the frame update fills in), plus the prediction group. The static
        reference is drawn from its own data and is not returned.
        """
        # the animated groups carry no data yet; the reference is drawn once and
        # goes last, on top, as it does in the static figures
        series = [(None, style) for style in roles]
        if ref_data is not None:
            series.append((ref_data, REFERENCE_STYLE))
        series.append((None, PREDICTION_STYLE))

        t_ref = t[0] + dt * np.arange(len(ref_data)) if ref_data is not None else None
        groups = _draw_panel(ax, t_ref, labels, ylabel, series, bounds=bounds, step=step)

        # frozen for the whole run, so the sliding window does not make the axes jump
        ax.set_ylim(*_ylim(ylim_data, pred_data, ref_data, np.asarray(bounds)))

        # hand back only the groups the frame update fills in
        return [group for group, (data, _) in zip(groups, series) if data is None]

    position_past, position_pred = setup_panel(
        ax_position, ['x', 'y', 'z'], 'position [m]',
        states['r'], pred_positions, ref[:, 0:3])

    velocity_past, velocity_pred = setup_panel(
        ax_velocity, ['u', 'v', 'w'], 'velocity (body) [m/s]',
        v_body, pred_v_body, ref_v_body)

    # the wind panel compares the estimate against the true wind, so it animates
    # two series; the prediction holds the estimate over the horizon
    band = n_sigma * wind_std
    wind_true_lines, wind_est_lines, wind_pred = setup_panel(
        ax_wind, ['$w_x$', '$w_y$', '$w_z$'], 'wind est. [m/s]',
        np.concatenate([wind_NED, wind_estimate - band, wind_estimate + band]),
        roles=(TRUE_STYLE, ESTIMATE_STYLE))

    dthrust_past, dthrust_pred = setup_panel(
        ax_dthrust, ['left', 'right'], 'dthrust [N/s]',
        U[:n_frames, 0:2], pred_dthrust, ref_controls[:, 0:2],
        bounds=DTHRUST_BOUNDS, step=True)

    # name the two regions the divider separates
    for x, text in ((1/3, 'past'), (5/6, 'MPC horizon')):
        ax_position.text(x, 1.02, text, transform=ax_position.transAxes,
                         ha='center', va='bottom', fontsize=9, color='gray')

    # 'now' divider and the shaded prediction region. The shading spans the full
    # height, so it is a rectangle in data-x / axes-y coordinates.
    dividers, shadings = [], []
    for ax in panels:
        dividers.append(ax.axvline(t[0], color='gray', linewidth=1))
        shadings.append(ax.add_patch(Rectangle(
            (t[0], 0), horizon, 1, color='gray', alpha=0.08, linewidth=0, zorder=0,
            transform=blended_transform_factory(ax.transData, ax.transAxes))))

    _style_legend(fig, _series_roles(reference=True, estimate=True, n_sigma=n_sigma)
                  + [('MPC prediction', PREDICTION_STYLE), ('bound', BOUND_STYLE)])

    title = fig.suptitle(f't = {t_sim[0]:.2f} s')
    wind_bands = []

    # ===== Frame update =====
    def update(frame):
        nonlocal plane_artists, wind_artists, wind_bands

        t_now = t_sim[frame]
        past = slice(0, frame + 1)
        t_past = t_sim[past]
        t_pred = t_now + dt * np.arange(N_horizon + 1)
        title.set_text(f't = {t_now:.2f} s')

        # --- 3D view ---
        for artist in plane_artists:
            artist.remove()
        plane_artists = _draw_plane(ax3d, states['r'][frame], states['R'][frame], scale=scale)

        trail.set_data(positions_ENU[past, 0], positions_ENU[past, 1])
        trail.set_3d_properties(positions_ENU[past, 2])

        # the corner arrow follows the wind of this frame, not its average
        for artist in wind_artists:
            artist.remove()
        wind_artists = _draw_wind_arrow(ax3d, wind_NED[frame], center, max_range)

        # --- Position and body velocity ---
        for k in range(3):
            position_past[k].set_data(t_past, states['r'][past, k])
            velocity_past[k].set_data(t_past, v_body[past, k])
            position_pred[k].set_data(t_pred, pred_positions[frame, :, k])
            velocity_pred[k].set_data(t_pred, pred_v_body[frame, :, k])

        # --- Wind estimate, held constant over the horizon ---
        for artist in wind_bands:
            artist.remove()
        wind_bands = []

        for k in range(3):
            estimate, deviation = wind_estimate[past, k], n_sigma * wind_std[past, k]
            wind_true_lines[k].set_data(t_past, wind_NED[past, k])
            wind_est_lines[k].set_data(t_past, estimate)
            wind_pred[k].set_data([t_now, t_now + horizon], [estimate[-1]] * 2)

            wind_bands.append(ax_wind.fill_between(
                t_past, estimate - deviation, estimate + deviation,
                color=f'C{k}', **BAND_STYLE))

        # --- dthrust, piecewise constant ---
        t_past_step = np.append(t_past, t_now + dt)
        for k in range(2):
            applied = U[past, k]
            predicted = pred_dthrust[frame, k]
            dthrust_past[k].set_data(t_past_step, np.append(applied, applied[-1]))
            dthrust_pred[k].set_data(t_pred, np.append(predicted, predicted[-1]))

        # --- Sliding window: the prediction keeps the right third, so the past
        #     grows leftwards out of the divider until it fills its two thirds ---
        panels[0].set_xlim(t_now - 2 * horizon, t_now + horizon)

        for divider, shading in zip(dividers, shadings):
            divider.set_xdata([t_now, t_now])
            shading.set_x(t_now)

        return []

    frames = range(0, n_frames, frame_step)
    interval = 1000 * dt * frame_step / time_scaling
    anim = FuncAnimation(fig, update, frames=frames, init_func=lambda: [],
                         interval=interval, repeat=False, blit=False)

    if output_path is None:
        return anim, fig

    fps = min(max(time_scaling / (dt * frame_step), 1), 60)
    anim.save(output_path, writer='pillow' if output_path.endswith('.gif') else 'ffmpeg',
              fps=fps, dpi=dpi)
    plt.close(fig)
    print(f"Animation saved to {output_path} ({len(frames)} frames, {fps:.1f} fps)")

    return None, None
