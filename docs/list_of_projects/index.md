<!-- # Publications and Projects that feature `acados`. -->
# Other Projects that feature `acados`

## Software interfaced with `acados`

- [openpilot](https://github.com/commaai/openpilot/)
is an open source driver assistance system.
[It has over 150 supported car makes and models.](https://github.com/commaai/openpilot/blob/master/docs/CARS.md)
`acados` is used within openpilot for lateral and longitudinal MPC.
It uses the `Cython` wrapper to the `acados` OCP solver in its software stack.

- [Rockit (Rapid Optimal Control kit)](https://gitlab.kuleuven.be/meco-software/rockit)
is a software framework to quickly prototype optimal control problems.
Notably, the software allows free end-time problems and multi-stage optimal problems.
The software is currently focused on direct methods and relies heavily on `CasADi`.
`acados` is interfaced as a `Rockit` solver by building on top of the Python interface of `acados`.

- [TuneMPC - a Python package for economic tuning of nonlinear model predictive control (NMPC) problems.](https://github.com/jdeschut/tunempc/)

- [OpenOCL](https://github.com/OpenOCL/OpenOCL)
is an open-source Matlab toolbox for modeling and solving optimal control problems.
It can use `CasADi` with IPOPT as a solver.
It also provides a higher level interface to `acados`, which is based on the Matlab interface of `acados`.

- [bioptim - a Python library for optimal control in biomechanics.](https://github.com/pyomeca/bioptim)

## Papers featuring `acados`
### with embedded deployment
<!-- in collaboration with syscop -->
- [Least Conservative Linearized Constraint Formulation for Real-Time Motion Generation](https://cdn.syscop.de/publications/Carlos2020.pdf)

- [An Efficient Real-Time NMPC for Quadrotor Position Control under Communication Time-Delay](https://cdn.syscop.de/publications/Carlos2020a.pdf)

- [NMPC for Racing Using a Singularity-Free Path-Parametric Model with Obstacle Avoidance](https://cdn.syscop.de/publications/Kloeser2020.pdf)

- [Mobility-enhanced MPC for Legged Locomotion on Rough Terrain](https://arxiv.org/abs/2105.05998)
    - [Video to Mobility-enhanced MPC for Legged Locomotion on Rough Terrain](https://www.dropbox.com/sh/mkr4pftcug6jlo7/AABNqu1AsGED2WSR8IqvaiUla?dl=0)

- [Continuous Control Set Nonlinear Model Predictive Control of Reluctance Synchronous Machines - IEEE Transactions on Control System Technology -- Andrea Zanelli et al 2021](https://ieeexplore.ieee.org/document/9360312)

- [Steering Action-aware Adaptive Cruise Control for Teleoperated Driving](https://ieeexplore.ieee.org/document/9945081)
  - [with public code on Github which has been applied on an F1TENTH vehicle](https://github.com/TUMFTM/tod_control/tree/2b67e8411e2ba1c5ddeb879d564ed28a989aebce/tod_shared_control)

### other
- [Contraction Properties of the Advanced Step Real-Time Iteration for NMPC at the IFAC World Congress 2020](https://cdn.syscop.de/publications/Nurkanovic2020b.pdf)

- [Real-Time Nonlinear Model Predictive Control for Microgrid Operation at the American Control Conference 2020](https://cdn.syscop.de/publications/Nurkanovic2020a.pdf)

- [Optimization-based Primary and Secondary Control of Microgrids](https://www.researchgate.net/profile/Armin_Nurkanovic/publication/341622767_Optimization-based_Primary_and_Secondary_Control_of_Microgrids/links/5f10519a299bf1e548ba5e77/Optimization-based-Primary-and-Secondary-Control-of-Microgrids.pdf)

- [TuneMPC—A Tool for Economic Tuning ofTracking (N)MPC Problems](https://cdn.syscop.de/publications/DeSchutter2020.pdf)

<!-- external -->
- [Model predictive control of wind turbine fatigue via online rainflow-counting on stress history and prediction](https://iopscience.iop.org/article/10.1088/1742-6596/1618/2/022041/pdf)

- [Embedded Real-Time Nonlinear Model Predictive Control for the Thermal Torque Derating of an Electric Vehicle, Winkler et al, IFAC 2021](https://cdn.syscop.de/publications/Winkler2021.pdf)

