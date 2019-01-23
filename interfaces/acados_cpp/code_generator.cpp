
#include "acados_cpp/code_generator.hpp"

namespace acados
{

code_generator::code_generator(ocp_nlp *nlp)
    : nlp_{nlp}
{
}

void code_generator::generate_s_function(std::string name)
{
    std::ofstream file;
    file.open("_autogen/" + name + ".c");

    generate_s_function_header(file, name);

    generate_mdl_initialize_sizes(file);

    generate_mdl_initialize_sample_times(file);

    generate_mdl_start(file);

    generate_mdl_outputs(file);

    generate_mdl_terminate(file);

    generate_s_function_footer(file);

    std::ofstream makefile;
    makefile.open("_autogen/" + name + "_make.m");

    generate_s_function_makefile(makefile, name);

    generate_dspace_makefile(makefile, name);
}

void code_generator::generate_s_function_header(std::ostream& out, std::string s_function_name)
{
    out << "\n#define S_FUNCTION_NAME " + s_function_name + "\n";
    out << "#define S_FUNCTION_LEVEL 2\n";

    out << "\n#include \"simstruc.h\"\n";

    out << "\n#include \"acados_c/external_function_interface.h\"\n";
    out << "#include \"acados_c/ocp_nlp_interface.h\"\n";

    if (nlp_->plan_->nlp_cost[0] == LINEAR_LS)
        out << "#include \"acados/ocp_nlp/ocp_nlp_cost_ls.h\"\n";
    else if (nlp_->plan_->nlp_cost[0] == NONLINEAR_LS)
        out << "#include \"acados/ocp_nlp/ocp_nlp_cost_nls.h\"\n";

    out << "\n#include \"blasfeo/include/blasfeo_d_aux.h\"\n";

    out << "\n";
    for (auto& module : nlp_->module_)
        out << "#include \"" + module.second.name() + ".h\"\n";

    out << "\n#define NUM_STAGES " + std::to_string(nlp_->N) + "\n";
    out << "#define LEN_INTERVAL " + std::to_string(nlp_->nlp_->Ts[0]) + "\n";
    out << "#define NUM_STATES " + std::to_string(nlp_->dims_->nx[0]) + "\n";
    out << "#define NUM_CONTROLS " + std::to_string(nlp_->dims_->nu[0]) + "\n";
}

void code_generator::generate_mdl_initialize_sizes(std::ostream& out)
{
    out << "\nstatic void mdlInitializeSizes(SimStruct *S)\n";
    out << "{\n";

    out << "\tssSetNumSFcnParams(S, 0);\n";
    out << "\tif (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))\n";
    out << "\t{\n";
    out << "\t\treturn; /* Parameter mismatch will be reported by Simulink */\n";
    out << "\t}\n";

    out << "\tif (!ssSetNumInputPorts(S, 1)) return;\n";
    out << "\tssSetInputPortWidth(S, 0, NUM_STATES);  /* x0 */\n";

    out << "\n\tssSetInputPortDirectFeedThrough(S, 0, true);\n";
    out << "\tssSetInputPortRequiredContiguous(S, 0, true);\n";

    out << "\n\tif (!ssSetNumOutputPorts(S, 2)) return;\n";
    out << "\tssSetOutputPortWidth(S, 0, NUM_CONTROLS);\n";
    out << "\tssSetOutputPortWidth(S, 1, 1);\n";

    out << "\n\tssSetNumPWork(S, 5);\n";

    out << "\n\tssSetNumSampleTimes(S, 1);\n";
    out << "}\n";
}

void code_generator::generate_mdl_initialize_sample_times(std::ostream& out)
{
    out << "\nstatic void mdlInitializeSampleTimes(SimStruct *S)\n";
    out << "{\n";
    out << "\tssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);\n";
    out << "\tssSetOffsetTime(S, 0, 0.0);\n";
    out << "}\n";
}

static void generate_ls_cost_initialization(std::ostream& out)
{
    out << "\n\tocp_nlp_cost_ls_model *stage_cost;\n";

    out << "\tfor (i = 0; i <= NUM_STAGES; ++i)\n";
    out << "\t{\n";
    out << "\t\tstage_cost = (ocp_nlp_cost_ls_model *) nlp_in->cost[i];\n";

    out << "\n\t\tblasfeo_dgese(nu[i] + nx[i], ny[i], 0.0, &stage_cost->Cyt, 0, 0);\n";
    out << "\t\tfor (j = 0; j < nu[i]; ++j)\n";
    out << "\t\t\tBLASFEO_DMATEL(&stage_cost->Cyt, j, nx[i] + j) = 1.0;\n";
    out << "\t\tfor (j = 0; j < nx[i]; ++j)\n";
    out << "\t\t\tBLASFEO_DMATEL(&stage_cost->Cyt, nu[i] + j, j) = 1.0;\n";

    out << "\n\t\tblasfeo_dgese(ny[i], ny[i], 0.0, &stage_cost->W, 0, 0);\n";
    out << "\t\tfor (j = 0; j < nx[i]; ++j)\n";
    out << "\t\t\tBLASFEO_DMATEL(&stage_cost->W, j, j) = 1.0;\n";
    out << "\t\tfor (j = 0; j < nu[i]; ++j)\n";
    out << "\t\t\tBLASFEO_DMATEL(&stage_cost->W, nx[i] + j, nx[i] + j) = 1.0;\n";

    out << "\n\t\tblasfeo_dvecse(nx[i], 0.0, &stage_cost->y_ref, 0);\n";
    out << "\t\tblasfeo_dvecse(nu[i], 0.0, &stage_cost->y_ref, nx[i]);\n";
    out << "\t}\n";
}

void code_generator::generate_mdl_start(std::ostream& out)
{
    out <<"\n#define MDL_START\n";
    out <<"static void mdlStart(SimStruct *S)\n";
    out << "{\n";
    out << "\tint nx[NUM_STAGES + 1], nu[NUM_STAGES + 1], ny[NUM_STAGES + 1],\n";
    out << "\t\tnb[NUM_STAGES + 1], nbx[NUM_STAGES + 1], nbu[NUM_STAGES + 1], \n";
    out << "\t\tng[NUM_STAGES + 1], nh[NUM_STAGES + 1], ns[NUM_STAGES + 1], nq[NUM_STAGES + 1],\n";
    out << "\t\tnz[NUM_STAGES + 1];\n";

    out << "\tint i, j;\n";
    out << "\n\tfor (i = 0; i < NUM_STAGES; ++i)\n";
    out << "\t{\n";
    out << "\t\tnx[i] = NUM_STATES;\n";
    out << "\t\tnu[i] = NUM_CONTROLS;\n";
    out << "\t\tny[i] = NUM_STATES + NUM_CONTROLS;\n";
    out << "\t\tnbx[i] = 0;\n";
    out << "\t\tnbu[i] = 0;\n";
    out << "\t\tnb[i] = nbx[i] + nbu[i];\n";
    out << "\t\tng[i] = 0;\n";
    out << "\t\tnh[i] = 0;\n";
    out << "\t\tns[i] = 0;\n";
    out << "\t\tnq[i] = 0;\n";
    out << "\t\tnz[i] = 0;\n";
    out << "\t}\n";

    out << "\n\tnbx[0] = NUM_STATES;\n";
    out << "\tnb[0] = nbx[0] + nbu[0];\n";

    out << "\n\tnx[NUM_STAGES] = NUM_STATES;\n";
    out << "\tnu[NUM_STAGES] = 0;\n";
    out << "\tny[NUM_STAGES] = NUM_STATES;\n";
    out << "\tnbx[NUM_STAGES] = 0;\n";
    out << "\tnbu[NUM_STAGES] = 0;\n";
    out << "\tng[NUM_STAGES] = 0;\n";
    out << "\tnh[NUM_STAGES] = 0;\n";
    out << "\tns[NUM_STAGES] = 0;\n";
    out << "\tnq[NUM_STAGES] = 0;\n";
    out << "\tnz[NUM_STAGES] = 0;\n";

    out << "\n\tocp_nlp_plan *plan = ocp_nlp_plan_create(NUM_STAGES);\n";
    out << "\tplan->nlp_solver = " + std::to_string(nlp_->plan_->nlp_solver) + ";\n";

    out << "\n\tfor (i = 0; i < NUM_STAGES; i++)\n";
    out << "\t{\n";
    out << "\t\tplan->nlp_cost[i] = " + std::to_string(nlp_->plan_->nlp_cost[0]) + ";\n";
    out << "\t\tplan->nlp_dynamics[i] = " + std::to_string(nlp_->plan_->nlp_dynamics[0]) + ";\n";
    out << "\t\tplan->sim_solver_plan[i].sim_solver = " +
            std::to_string(nlp_->plan_->sim_solver_plan[0].sim_solver) + ";\n";
    out << "\t\tplan->nlp_constraints[i] = 0;\n";
    out << "\t}\n";

    out << "\tplan->nlp_constraints[NUM_STAGES] = 0;\n";
    out << "\tplan->nlp_cost[NUM_STAGES] = " + std::to_string(nlp_->plan_->nlp_cost[nlp_->N])
                                             + ";\n";

    out << "\tplan->ocp_qp_solver_plan.qp_solver = " +
            std::to_string(nlp_->plan_->ocp_qp_solver_plan.qp_solver) + ";\n";

    out << "\n\tocp_nlp_config *config = ocp_nlp_config_create(*plan);\n";

    out << "\n\tocp_nlp_dims *nlp_dims = ocp_nlp_dims_create(config);\n";
    // TODO(oj): this cant work anymore, check if can be fixed
    out << "\tocp_nlp_dims_initialize(config, nx, nu, ny, nbx, nbu, ng, nh, nq, ns, nz, ";
    out << "nlp_dims);\n";

    out << "\n\texternal_function_casadi *expl_vde_for = malloc(NUM_STAGES "
           "* sizeof(external_function_casadi));\n";
    out << "\tfor (i = 0; i < NUM_STAGES; ++i)\n";
    out << "\t{\n";
    out << "\t\texpl_vde_for[i].casadi_fun = &" + nlp_->cached_model_ + ";\n";
    out << "\t\texpl_vde_for[i].casadi_n_in = &" + nlp_->cached_model_ + "_n_in;\n";
    out << "\t\texpl_vde_for[i].casadi_n_out = &" + nlp_->cached_model_ + "_n_out;\n";
    out << "\t\texpl_vde_for[i].casadi_sparsity_in = &" + nlp_->cached_model_ + "_sparsity_in;\n";
    out << "\t\texpl_vde_for[i].casadi_sparsity_out = &" + nlp_->cached_model_ + "_sparsity_out;\n";
    out << "\t\texpl_vde_for[i].casadi_work = &" + nlp_->cached_model_ + "_work;\n";
    out << "\t}\n";
    out << "\texternal_function_casadi_create_array(NUM_STAGES, expl_vde_for);\n";

    out << "\n\tocp_nlp_in *nlp_in = ocp_nlp_in_create(config, nlp_dims);\n";

    out << "\n\tfor (i = 0; i < NUM_STAGES; ++i)\n";
    out << "\t\tnlp_in->Ts[i] = LEN_INTERVAL;\n";

    out << "\n\tfor (i = 0; i < NUM_STAGES; ++i)\n";
    out << "\t\tocp_nlp_dynamics_model_set(config, nlp_in, i, \"expl_vde_for\"";
    out << ", &expl_vde_for[i]);\n";

    generate_ls_cost_initialization(out);

    out << "\n\tocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) "
           "nlp_in->constraints;\n";

    out << "\tfor (i = 0; i < NUM_STATES; ++i)\n";
    out << "\t\tconstraints[0]->idxb[i] = NUM_CONTROLS + i;\n";

    out << "\n\tvoid *nlp_opts = ocp_nlp_opts_create(config, nlp_dims);\n";

    out << "\n\tocp_nlp_out *nlp_out = ocp_nlp_out_create(config, nlp_dims);\n";

    out << "\n\tocp_nlp_solver *nlp_solver = ocp_nlp_solver_create(config, nlp_dims, nlp_opts);\n";

    out << "\n\tssGetPWork(S)[0] = (void *) nlp_dims;\n";
    out << "\tssGetPWork(S)[1] = (void *) nlp_in;\n";
    out << "\tssGetPWork(S)[2] = (void *) nlp_out;\n";
    out << "\tssGetPWork(S)[3] = (void *) nlp_opts;\n";
    out << "\tssGetPWork(S)[4] = (void *) nlp_solver;\n";
    out << "}\n";
}

void code_generator::generate_mdl_outputs(std::ostream& out)
{
    out << "\nstatic void mdlOutputs(SimStruct *S, int_T tid)\n";
    out << "{\n";
    out << "\tocp_nlp_in *nlp_in = (ocp_nlp_in *) ssGetPWork(S)[1];\n";
    out << "\tocp_nlp_out *nlp_out = (ocp_nlp_out *) ssGetPWork(S)[2];\n";
    out << "\tocp_nlp_solver *nlp_solver = (ocp_nlp_solver *) ssGetPWork(S)[4];\n";

    out << "\n\tocp_nlp_constraints_bgh_model **constraints = (ocp_nlp_constraints_bgh_model **) "
           "nlp_in->constraints;\n";

    out << "\n\tconst double *x0 = ssGetInputPortRealSignal(S, 0);\n";
    out << "\tblasfeo_pack_dvec(NUM_STATES, (double *) x0, &constraints[0]->d, 0);\n";
    out << "\tblasfeo_pack_dvec(NUM_STATES, (double *) x0, &constraints[0]->d, NUM_STATES);\n";

    out << "\n\tint status = ocp_nlp_solve(nlp_solver, nlp_in, nlp_out);\n";

    out << "\n\tdouble *u0_opt = ssGetOutputPortRealSignal(S, 0);\n";
    out << "\tdouble *status_out = ssGetOutputPortRealSignal(S, 1);\n";

    out << "\n\tblasfeo_unpack_dvec(NUM_CONTROLS, &nlp_out->ux[0], 0, u0_opt);\n";
    out << "\t*status_out = (double) status;\n";
    out << "}\n";
}

void code_generator::generate_mdl_terminate(std::ostream& out)
{
    out << "\nstatic void mdlTerminate(SimStruct *S)\n";
    out << "{\n";
    out << "\tfree(ssGetPWork(S)[0]);\n";
    out << "\tfree(ssGetPWork(S)[1]);\n";
    out << "\tfree(ssGetPWork(S)[2]);\n";
    out << "\tfree(ssGetPWork(S)[3]);\n";
    out << "\tfree(ssGetPWork(S)[4]);\n";
    out << "}\n";
}

void code_generator::generate_s_function_footer(std::ostream& out)
{
    out << "\n#ifdef MATLAB_MEX_FILE /* Is this file being compiled as a MEX-file? */\n";
    out << "#include \"simulink.c\"  /* MEX-file interface mechanism */\n";
    out << "#else\n";
    out << "#include \"cg_sfun.h\" /* Code generation registration function */\n";
    out << "#endif\n";
}

void code_generator::generate_s_function_makefile(std::ostream& out, std::string function_name)
{
    out << "\n";

    out << "% Dialog with which the user selects the folder where the acados libs\n";
    out << "% reside.\n";
    out << "\nif(~ exist('acados_lib_path', 'var'))\n";
    out << "\tacados_path = uigetdir('', 'Please select folder with the acados libraries.');";
    out << "\n";
    out << "\tacados_lib_path = fullfile(acados_path, 'lib');\n";
    out << "\tacados_include_path = fullfile(acados_path, 'include');\n";
    out << "\tblasfeo_include_path = fullfile(acados_path, 'include/blasfeo/include');\n";
    out << "end\n";

    out << "\nmex_command = 'mex " + function_name + ".c';\n";
    out << "mex_command = [mex_command, ' " + nlp_->cached_model_ + ".c'];\n";
    out << "mex_command = [mex_command, ' -I', acados_include_path];\n";
    out << "mex_command = [mex_command, ' -I', blasfeo_include_path];\n";
    out << "mex_command = [mex_command, ' -L', acados_lib_path];\n";
    out << "mex_command = [mex_command, ' -lacados -lhpmpc -lhpipm -lqpOASES_e -lblasfeo'];\n";

    out << "\neval(mex_command);\n";
}

void code_generator::generate_dspace_makefile(std::ostream& out, std::string function_name)
{
    out << "\n";

    out << "model_name = 'acados_template';\n";

    out << "if isempty(which(model_name))\n";
    out << "\tnew_system(model_name);\n";
    out << "else\n";
    out << "\tload_system(model_name);\n";
    out << "end\n";

    out << "\nsave_system(model_name);\n";

    out << "\nif (isempty(find_system(model_name, 'Name', 'S-function')))\n";
    out << "\tadd_block('simulink/User-Defined Functions/S-Function', "
           "[model_name, '/S-function']);\n";
    out << "\tset_param([model_name, '/S-function'], 'FunctionName', '" + function_name + "');\n";
    out << "\tset_param([model_name, '/S-function'], 'position', [100 40 220 100]);\n";
    out << "end\n";

    out << "\nif (isempty(find_system(model_name, 'Name', 'u_opt')))\n";
    out << "\tadd_block('simulink/Sinks/Scope', [model_name, '/u_opt']);\n";
    out << "\tset_param([model_name, '/u_opt'], 'position', [250 40 270 60]);\n";
    out << "\tadd_line(model_name, 'S-function/1', 'u_opt/1');\n";
    out << "end\n";

    out << "\nif (isempty(find_system(model_name, 'Name', 'Status')))";
    out << "\tadd_block('simulink/Sinks/Scope', [model_name, '/Status']);\n";
    out << "\tset_param([model_name, '/Status'], 'position', [250 80 270 100]);\n";
    out << "\tadd_line(model_name, 'S-function/2', 'Status/1');\n";
    out << "end\n";

    out << "\nset_param(getActiveConfigSet(model_name), 'StopTime', 'Inf');\n";
    out << "\nset_param(getActiveConfigSet(model_name), 'Solver', 'FixedStepDiscrete');\n";
    out << "\nset_param(getActiveConfigSet(model_name), 'FixedStep', '0.1');\n";

    out << "\nset_param(getActiveConfigSet(model_name), 'CustomSource', [pwd, '/" +
           nlp_->cached_model_ + ".c']);\n";
    out << "set_param(getActiveConfigSet(model_name), 'CustomInclude', "
           "[acados_include_path, ' ', blasfeo_include_path]);\n";
    out << "set_param(getActiveConfigSet(model_name), 'CustomLibrary', ...\n";
    out << "\t[acados_lib_path, '/../dspace/lib/acados.lib ', ...\n";
    out << "\tacados_lib_path, '/../dspace/lib/hpipm.lib ', ...\n";
    out << "\tacados_lib_path, '/../dspace/lib/hpmpc.lib ', ...\n";
    out << "\tacados_lib_path, '/../dspace/lib/blasfeo.lib ']);\n";

    out << "\nsave_system(model_name);\n";

    out << "\nopen_system(model_name);\n";

}

}  // namespace acados
