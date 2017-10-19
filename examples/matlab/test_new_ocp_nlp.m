import acados.*
[ode_fun, nx, nu] = chen_model();

nlpfun = ocp_nlp_function(ode_fun);