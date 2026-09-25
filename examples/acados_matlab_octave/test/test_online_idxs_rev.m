import casadi.*

N = 20;
T_HORIZON = 1.0;
X0 = [2.0; 0.0];
V_MAX = 1.5;
U_MAX = 3.0;

model = AcadosModel();
p = SX.sym('p');
v = SX.sym('v');
u = SX.sym('u');
model.name = 'online_idxs_rev';
model.x = vertcat(p, v);
model.u = u;
model.xdot = SX.sym('xdot', 2);
model.f_expl_expr = vertcat(v, u);
model.f_impl_expr = model.f_expl_expr - model.xdot;

ocp = AcadosOcp();
ocp.model = model;
ocp.solver_options.N_horizon = N;
ocp.solver_options.tf = T_HORIZON;
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.nlp_solver_type = 'SQP';
ocp.solver_options.nlp_solver_max_iter = 1;
ocp.solver_options.nlp_solver_tol_stat = 1e-10;
ocp.solver_options.nlp_solver_tol_eq = 1e-10;
ocp.solver_options.nlp_solver_tol_ineq = 1e-10;
ocp.solver_options.nlp_solver_tol_comp = 1e-10;
ocp.solver_options.eval_residual_at_max_iter = true;
ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM';

ocp.cost.cost_type = 'NONLINEAR_LS';
ocp.cost.W = diag([1.0, 1e-1, 1e-3]);
ocp.cost.yref = zeros(3, 1);
ocp.model.cost_y_expr = vertcat(model.x, model.u);

ocp.constraints.lbu = -U_MAX;
ocp.constraints.ubu = U_MAX;
ocp.constraints.idxbu = 0;
ocp.constraints.idxbx = 1;
ocp.constraints.lbx = -V_MAX;
ocp.constraints.ubx = V_MAX;
ocp.constraints.x0 = X0;

ocp.constraints.idxs_rev_0 = [0; -1; -1];
ocp.constraints.idxs_rev = [0; 1];
ocp.cost.zl_0 = 0;
ocp.cost.Zl_0 = 2;
ocp.cost.zu_0 = 0;
ocp.cost.Zu_0 = 2;
ocp.cost.zl = zeros(2, 1);
ocp.cost.Zl = 2 * ones(2, 1);
ocp.cost.zu = zeros(2, 1);
ocp.cost.Zu = 2 * ones(2, 1);

ocp_solver = AcadosOcpSolver(ocp);
variants = {'individual', 'all hard', 'joint', 'u hard', 'x hard'};
unused_slacks = {[], [1, 2], 2, 2, 2};

for variant_idx = 1:length(variants)
    variant = variants{variant_idx};

    switch variant
        case 'individual'
            idxs_rev_0 = [0; -1; -1];
            idxs_rev = [0; 1];
        case 'all hard'
            idxs_rev_0 = [-1; -1; -1];
            idxs_rev = [-1; -1];
        case 'joint'
            idxs_rev_0 = [0; -1; -1];
            idxs_rev = [0; 0];
        case 'u hard'
            idxs_rev_0 = [-1; -1; -1];
            idxs_rev = [-1; 0];
        case 'x hard'
            idxs_rev_0 = [0; -1; -1];
            idxs_rev = [0; -1];
    end

    ocp_solver.set('constr_idxs_rev', idxs_rev_0, 0);
    for stage = 1:N-1
        ocp_solver.set('constr_idxs_rev', idxs_rev, stage);
    end

    ocp_solver.solve();
    ocp_solver.print_statistics();
    status = ocp_solver.get('status');
    assert(status == 0, sprintf('variant %s returned status %d', variant, status));

    input_traj = zeros(N, 1);
    for stage = 0:N-1
        input_traj(stage+1) = ocp_solver.get('u', stage);
    end
    if any(strcmp(variant, {'all hard', 'u hard'}))
        assert(all(input_traj < U_MAX), 'hard input upper bound violated');
        assert(all(-U_MAX < input_traj), 'hard input lower bound violated');
    else
        assert(any(input_traj >= U_MAX) && any(input_traj <= -U_MAX), ...
            'input bound softening should be exploited, otherwise test problem is bad');
    end

    velocity_traj = zeros(N-1, 1);
    for stage = 1:N-1
        state = ocp_solver.get('x', stage);
        velocity_traj(stage) = state(2);
    end
    if any(strcmp(variant, {'all hard', 'x hard'}))
        assert(all(velocity_traj < V_MAX), 'hard state upper bound violated');
        assert(all(-V_MAX < velocity_traj), 'hard state lower bound violated');
    else
        assert(any(velocity_traj >= V_MAX) || any(velocity_traj <= -V_MAX), ...
            'state bound softening should be exploited, otherwise test problem is bad');
    end

    nb = ocp.dims.nbx + ocp.dims.nbu;
    ns = ocp.dims.ns;
    for stage = 1:N-1
        sl = ocp_solver.get('sl', stage);
        su = ocp_solver.get('su', stage);
        lam = ocp_solver.get('lam', stage);
        lam_l = lam(1:nb);
        lam_u = lam(nb+1:2*nb);
        lam_sl = lam(2*nb+1:2*nb+ns);
        lam_su = lam(2*nb+ns+1:2*nb+2*ns);
        for slack_idx = unused_slacks{variant_idx}
            assert(abs(sl(slack_idx)) < 1e-8, ...
                sprintf('unused sl(%d) is nonzero for %s', slack_idx, variant));
            assert(abs(su(slack_idx)) < 1e-8, ...
                sprintf('unused su(%d) is nonzero for %s', slack_idx, variant));
            assert(abs(lam_sl(slack_idx)) < 1e-8, ...
                sprintf('unused lower slack multiplier is nonzero for %s at stage %d', variant, stage));
            assert(abs(lam_su(slack_idx)) < 1e-8, ...
                sprintf('unused upper slack multiplier is nonzero for %s at stage %d', variant, stage));
        end

        if strcmp(variant, 'joint')
            input_value = ocp_solver.get('u', stage);
            state = ocp_solver.get('x', stage);
            velocity_value = state(2);

            max_lower_violation = max([-U_MAX-input_value, -V_MAX-velocity_value, 0.0]);
            max_upper_violation = max([input_value-U_MAX, velocity_value-V_MAX, 0.0]);
            tol_slack = 5e-5;
            lower_slack_error = abs(sl(1) - max_lower_violation);
            upper_slack_error = abs(su(1) - max_upper_violation);
            assert(lower_slack_error <= tol_slack * (1 + abs(max_lower_violation)), ...
                sprintf('lower joint slack mismatch at stage %d', stage));
            assert(upper_slack_error <= tol_slack * (1 + abs(max_upper_violation)), ...
                sprintf('upper joint slack mismatch at stage %d', stage));

            cost_scale = ocp.solver_options.cost_scaling(stage+1);
            sl_stat = cost_scale * (ocp.cost.zl(1) + ocp.cost.Zl(1)*sl(1)) ...
                - sum(lam_l) - lam_sl(1);
            su_stat = cost_scale * (ocp.cost.zu(1) + ocp.cost.Zu(1)*su(1)) ...
                - sum(lam_u) - lam_su(1);
            assert(abs(sl_stat) <= 1e-10 && abs(su_stat) <= 1e-10, ...
                sprintf('slack stationarity failed at stage %d: lower=%g upper=%g', ...
                    stage, sl_stat, su_stat));
        end
    end
end

clear ocp_solver
disp('online idxs_rev updates passed');
