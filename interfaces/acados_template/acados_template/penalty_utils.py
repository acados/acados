from typing import Optional

from .acados_ocp import AcadosOcp
from scipy.linalg import block_diag

import casadi as ca
import numpy as np


def formulate_constraints_as_penalty(
    ocp: AcadosOcp,
    constr_expr: ca.SX,
    weight: float,
    upper_bound: Optional[float],
    lower_bound: Optional[float],
    penalty_type: Optional[str] = "L2",
) -> AcadosOcp:

    if upper_bound is None and lower_bound is None:
        raise ValueError("Either upper or lower bound must be provided.")

    violation_expr = 0.0
    if upper_bound is not None:
        violation_expr = ca.fmax(violation_expr, (constr_expr - upper_bound))
    if lower_bound is not None:
        violation_expr = ca.fmax(violation_expr, (lower_bound - constr_expr))

    if penalty_type == "L2":
        ocp.cost.yref = np.concatenate((ocp.cost.yref, np.zeros(1)))
        ocp.model.cost_y_expr = ca.vertcat(ocp.model.cost_y_expr, violation_expr)
        if ocp.cost.cost_type == "NONLINEAR_LS":
            ocp.cost.W = block_diag(ocp.cost.W, weight)
        elif ocp.cost.cost_type == "CONVEX_OVER_NONLINEAR":
            new_residual = ca.SX.sym('new_residual', violation_expr.shape)
            ocp.model.cost_r_in_psi_expr = ca.vertcat(ocp.model.cost_r_in_psi_expr, new_residual)
            ocp.model.cost_psi_expr += .5 * weight * new_residual**2
        else:
            raise NotImplementedError(
                f"Penalty type {penalty_type} is not yet supported for ocp.."
            )
    else:
        raise NotImplementedError(f"Penalty type {penalty_type} is not yet supported.")

    return ocp
