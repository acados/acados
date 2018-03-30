# Author: Gianluca Frison



include ../../Makefile.rule



LIBS =

ifeq ($(ACADOS_WITH_C_INTERFACE), 1)
LIBS += -lacados_c
endif

LIBS += -lacore -lhpipm

ifeq ($(ACADOS_WITH_QPOASES), 1)
LIBS += -lqpOASES_e
endif

ifeq ($(ACADOS_WITH_HPMPC), 1)
LIBS += -lhpmpc
endif

ifeq ($(ACADOS_WITH_QPDUNES), 1)
CCC = $(CXX)
LIBS += -lqpdunes -lqpOASES
else
CCC = $(CC)
endif

ifeq ($(ACADOS_WITH_QORE), 1)
LIBS += -lqore_dense
endif

LIBS += -lblasfeo -lm -lblas -llapack



# Comment this out to enable using gprof
# CFLAGS  += -pg
# LDFLAGS += -pg

EXAMPLES =
EXAMPLES += sim_wt_model
EXAMPLES += sim_wt_model_nx6
EXAMPLES += sim_crane_example
EXAMPLES += mass_spring_example
#EXAMPLES += mass_spring_pcond_hpipm_split
#EXAMPLES += mass_spring_fcond_hpipm_split
#EXAMPLES += mass_spring_offline_fcond_qpoases_split
EXAMPLES += nonlinear_chain_ocp_nlp
EXAMPLES += sim_crane_no_interface
EXAMPLES += mass_spring_example_no_interface
EXAMPLES += nonlinear_chain_ocp_nlp_no_interface
EXAMPLES += wind_turbine_nmpc

examples: $(EXAMPLES)



RUN_EXAMPLES =
RUN_EXAMPLES += run_sim_wt_model
RUN_EXAMPLES += run_sim_wt_model_nx6
RUN_EXAMPLES += run_sim_crane_example
RUN_EXAMPLES += run_mass_spring_example
#RUN_EXAMPLES += run_mass_spring_pcond_hpipm_split
#RUN_EXAMPLES += run_mass_spring_fcond_hpipm_split
#RUN_EXAMPLES += run_mass_spring_offline_fcond_qpoases_split
RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp
RUN_EXAMPLES += run_sim_crane_no_interface
RUN_EXAMPLES += run_mass_spring_example_no_interface
RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp_no_interface
RUN_EXAMPLES += run_wind_turbine_nmpc

run_examples: $(RUN_EXAMPLES)



#################################################
# crane model
#################################################

CRANE_OBJS =
CRANE_OBJS += crane_model/ode_model.o
CRANE_OBJS += crane_model/jac_model.o
CRANE_OBJS += crane_model/vde_adj_model.o
CRANE_OBJS += crane_model/vde_forw_model.o
CRANE_OBJS += crane_model/vde_hess_model.o
CRANE_OBJS += crane_model/impl_ode.o
CRANE_OBJS += crane_model/impl_ode_jac_x.o
CRANE_OBJS += crane_model/impl_ode_jac_xdot.o
CRANE_OBJS += crane_model/impl_ode_jac_u.o
CRANE_OBJS += crane_model/impl_ode_inc_J_x_xdot.o
CRANE_OBJS += crane_model/impl_ode_J_x_u.o
CRANE_OBJS += crane_model/impl_ode_J_x_xdot_u.o

sim_crane_no_interface: $(CRANE_OBJS) sim_crane_no_interface.o
	$(CCC) -o sim_crane_no_interface.out sim_crane_no_interface.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_crane_no_interface build complete."
	@echo

run_sim_crane_no_interface:
	./sim_crane_no_interface.out



sim_crane_example: $(CRANE_OBJS) sim_crane_example.o
	$(CCC) -o sim_crane_example.out sim_crane_example.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_crane_example build complete."
	@echo

run_sim_crane_example:
	./sim_crane_example.out



#################################################
# wind turbine model with 6 states
#################################################

WT_NX6_OBJS =
WT_NX6_OBJS += wt_model_nx6/expl_ode_fun.o
WT_NX6_OBJS += wt_model_nx6/expl_vde_for.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_fun.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_x.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_xdot.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_u.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_fun_jac_x_xdot.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_x_xdot_u.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_x_u.o

sim_wt_model_nx6: $(WT_NX6_OBJS) sim_wt_model_nx6.o
	$(CCC) -o sim_wt_model_nx6.out sim_wt_model_nx6.o  $(WT_NX6_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_wt_model_nx6 build complete."
	@echo

run_sim_wt_model_nx6:
	./sim_wt_model_nx6.out



WT_NX6P2_OBJS =
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/expl_ode_fun.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/expl_vde_for.o

wind_turbine_nmpc: $(WT_NX6P2_OBJS) wind_turbine_nmpc.o
	$(CCC) -o wind_turbine_nmpc.out wind_turbine_nmpc.o  $(WT_NX6P2_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example wind_turbine_nmpc build complete."
	@echo

run_wind_turbine_nmpc:
	./wind_turbine_nmpc.out



#################################################
# wind turbine model with 3 states
#################################################

WT_OBJS =
WT_OBJS += wt_model_nx3/ode_energy_balanced_model.o
WT_OBJS += wt_model_nx3/jac_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_adj_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_hess_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_odeFun_energy_balanced_model.o
#WT_OBJS += wt_model_nx3/impl_jacFun_x_energy_balanced_model.o
#WT_OBJS += wt_model_nx3/impl_jacFun_xdot_energy_balanced_model.o
#WT_OBJS += wt_model_nx3/impl_jacFun_u_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_ode_inc_J_x_xdot_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_ode_J_x_u_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_ode_J_x_xdot_u_energy_balanced_model.o

sim_wt_model: $(WT_OBJS) sim_wt_model.o
	$(CCC) -o sim_wt_model.out sim_wt_model.o  $(WT_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_wt_model build complete."
	@echo

run_sim_wt_model:
	./sim_wt_model.out



#################################################
# (linear) mass spring model
#################################################

MASS_SPRING_OBJS =
MASS_SPRING_OBJS += mass_spring_model/mass_spring_qp.o

mass_spring_example_no_interface: $(MASS_SPRING_OBJS) mass_spring_example_no_interface.o
	$(CCC) -o mass_spring_example_no_interface.out $(MASS_SPRING_OBJS) mass_spring_example_no_interface.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example_no_interface build complete."
	@echo

run_mass_spring_example_no_interface:
	./mass_spring_example_no_interface.out



mass_spring_example: $(MASS_SPRING_OBJS) mass_spring_example.o
	$(CCC) -o mass_spring_example.out $(MASS_SPRING_OBJS) mass_spring_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example build complete."
	@echo

run_mass_spring_example:
	./mass_spring_example.out



mass_spring_pcond_hpipm_split: mass_spring_pcond_hpipm_split.o
	$(CCC) -o mass_spring_pcond_hpipm_split.out mass_spring_pcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_hpipm_split build complete."
	@echo

run_mass_spring_pcond_hpipm_split:
	./mass_spring_pcond_hpipm_split.out



mass_spring_fcond_hpipm_split: mass_spring_fcond_hpipm_split.o
	$(CCC) -o mass_spring_fcond_hpipm_split.out mass_spring_fcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_fcond_hpipm_split build complete."
	@echo

run_mass_spring_fcond_hpipm_split:
	./mass_spring_fcond_hpipm_split.out



mass_spring_offline_fcond_qpoases_split: mass_spring_offline_fcond_qpoases_split.o
	$(CCC) -o mass_spring_offline_fcond_qpoases_split.out mass_spring_offline_fcond_qpoases_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_offline_fcond_qpoases_split build complete."
	@echo

run_mass_spring_offline_fcond_qpoases_split:
	./mass_spring_offline_fcond_qpoases_split.out



#################################################
# (nonlinear) chain model
#################################################

CHAIN_OBJS =
### explicit ###
# 2 masses
CHAIN_OBJS += chain_model/vde_chain_nm2.o
CHAIN_OBJS += chain_model/jac_chain_nm2.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm2.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm2.o
CHAIN_OBJS += chain_model/ls_cost_nm2.o
CHAIN_OBJS += chain_model/ls_costN_nm2.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm2.o
# 3 masses
CHAIN_OBJS += chain_model/vde_chain_nm3.o
CHAIN_OBJS += chain_model/jac_chain_nm3.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm3.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm3.o
CHAIN_OBJS += chain_model/ls_cost_nm3.o
CHAIN_OBJS += chain_model/ls_costN_nm3.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm3.o
# 4 masses
CHAIN_OBJS += chain_model/vde_chain_nm4.o
CHAIN_OBJS += chain_model/jac_chain_nm4.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm4.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm4.o
CHAIN_OBJS += chain_model/ls_cost_nm4.o
CHAIN_OBJS += chain_model/ls_costN_nm4.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm4.o
# 5 masses
CHAIN_OBJS += chain_model/vde_chain_nm5.o
CHAIN_OBJS += chain_model/jac_chain_nm5.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm5.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm5.o
CHAIN_OBJS += chain_model/ls_cost_nm5.o
CHAIN_OBJS += chain_model/ls_costN_nm5.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm5.o
# 6 masses
CHAIN_OBJS += chain_model/vde_chain_nm6.o
CHAIN_OBJS += chain_model/jac_chain_nm6.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm6.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm6.o
CHAIN_OBJS += chain_model/ls_cost_nm6.o
CHAIN_OBJS += chain_model/ls_costN_nm6.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm6.o
### implicit ###
# 2 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_xdot_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_u_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_u_chain_nm2.o
# 3 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_xdot_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_u_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_u_chain_nm3.o
# 4 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_xdot_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_u_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_u_chain_nm4.o
# 5 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_xdot_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_u_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_u_chain_nm5.o
# 6 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_xdot_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_u_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_u_chain_nm6.o

nonlinear_chain_ocp_nlp_no_interface: $(CHAIN_OBJS) nonlinear_chain_ocp_nlp_no_interface.o
	$(CCC) -o nonlinear_chain_ocp_nlp_no_interface.out $(CHAIN_OBJS) nonlinear_chain_ocp_nlp_no_interface.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp_no_interface build complete."
	@echo

run_nonlinear_chain_ocp_nlp_no_interface:
	./nonlinear_chain_ocp_nlp_no_interface.out



nonlinear_chain_ocp_nlp: $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o
	$(CCC) -o nonlinear_chain_ocp_nlp.out $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo

run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out



#################################################
# clean
#################################################

clean:
	rm -f *.o
	rm -f *.out

deep_clean: clean
	rm -f chain_model/*.o
	rm -f wt_model_nx3/*.o
	rm -f wt_model_nx6/*.o
	rm -f wt_model_nx6/nx6p2/*.o

