# Periodic trajectory tracking for a fixed-wing aircraft

``` eval_rst
.. meta::
   :description: Periodic trajectory tracking for a fixed-wing aircraft using nonlinear MPC and an EKF for state and wind estimation
   :keywords: fixed-wing aircraft, periodic trajectory tracking, wind estimation, nonlinear model predictive control, documentation example
```

This example uses nonlinear MPC to make a fixed-wing aircraft track a periodic reference, a tilted circular orbit, in unknown and time-varying wind.
An extended Kalman filter (EKF) estimates the state and the wind.
The control scheme follows a very similar approach to [Gros et al. (2012)](https://publications.syscop.de/Gros2012e.pdf).
In [Harzer et al. (2026)](https://www.imavs.org/pdf/imav.2026.6), the same approach is applied to a tethered aircraft, combined with a more complex wind field estimation.

The example is in [`examples/acados_python/fixed_wing_aircraft`](https://github.com/acados/acados/tree/main/examples/acados_python/fixed_wing_aircraft).

``` eval_rst
.. raw:: html

   <video src="../../_static/fixed_wing_aircraft_example_animation_3d_with_states.mp4" width="100%" autoplay loop muted playsinline controls
          title="Closed-loop simulation of the fixed-wing aircraft tracking a tilted circular orbit"></video>
```

## Model

The aircraft is modeled as a rigid body, the aerodynamic forces and moments come from a polynomial regression model.
Its 22-dimensional state consists of position (3 dimensions), velocity (3), orientation as a rotation matrix (9), rotation rates (3), left and right thrust (2), and elevator and aileron deflection (2).
The 4 controls are the rates of the left thrust, the right thrust, the elevator deflection and the aileron deflection.
While the two thrusts act independently, the ailerons deflect anti-symmetrically, so one deflection describes both.
Wind enters as a 3-dimensional parameter and is assumed constant over the prediction horizon.

## Optimal control problem

```math
\begin{aligned}
&\underset{\substack{x_0, \dots, x_N \\ u_0, \dots, u_{N-1}}}{\mathrm{minimize}}
&& \frac{1}{2} \sum_{k=0}^{N} \left\| x_k - x_k^\mathrm{ref} \right\|_Q^2
+ \frac{1}{2} \sum_{k=0}^{N-1} \left\| u_k - u_k^\mathrm{ref} \right\|_R^2 \\[4ex]
&\text{subject to}
&& \begin{aligned}[t]
0 &= F(x_{k+1}, x_k, u_k, p), && k = 0, \dots, N-1, \\[1ex]
x_0 &= \hat x, \\[1ex]
\underline{u} \le u_k &\le \overline{u}, && k = 0, \dots, N-1, \\[1ex]
\underline{x} \le x_k &\le \overline{x}, && k = 0, \dots, N, \\[1ex]
\underline{h} \le h(x_k) &\le \overline{h}, && k = 0, \dots, N,
\end{aligned}
\end{aligned}
```

where `$x_k^\mathrm{ref}$` and `$u_k^\mathrm{ref}$` are a precalculated reference trajectory, `$\hat x$` the estimated current state and `$p$` the estimated wind.
`$h$` are nonlinear constraints consisting of the airspeed, the angle of attack and the sideslip angle. They are optionally softened by slacks.

## State and wind estimation

The EKF estimates the full state and, in addition, the 3-dimensional wind parameter.
The prediction step integrates the dynamics with an `AcadosSimSolver`.
Its forward sensitivities with respect to the state and the parameter yield the required Jacobian.
