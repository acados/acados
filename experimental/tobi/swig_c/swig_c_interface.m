import acados.*

plan = ocp_nlp_plan_create(NN);

plan.nlp_solver = SQP;

    for i=1:NN
        plan->nlp_cost[i] = LINEAR_LS;

    plan->ocp_qp_solver_plan.qp_solver = PARTIAL_CONDENSING_HPIPM;


    for (int i = 0; i < NN; i++)
    {
        plan->nlp_dynamics[i] = CONTINUOUS_MODEL;

        plan->sim_solver_plan[i].sim_solver = GNSF;
    }

    for (int i = 0; i <= NN; i++)
        plan->nlp_constraints[i] = BGH;

ocp_nlp_solver_config *config = ocp_nlp_config_create(*plan);