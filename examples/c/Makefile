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
EXAMPLES += sim_gnsf_crane
EXAMPLES += sim_gnsf_crane_empty_z
EXAMPLES += sim_irk_crane2
EXAMPLES += sim_wt_model
EXAMPLES += sim_wt_model_nx6
EXAMPLES += sim_crane_example
EXAMPLES += mass_spring_example
#EXAMPLES += mass_spring_pcond_hpipm_split
#EXAMPLES += mass_spring_fcond_hpipm_split
#EXAMPLES += mass_spring_offline_fcond_qpoases_split
EXAMPLES += nonlinear_chain_ocp_nlp
EXAMPLES += sim_crane_no_interface
# EXAMPLES += mass_spring_example_no_interface
EXAMPLES += nonlinear_chain_ocp_nlp_no_interface
EXAMPLES += sim_gnsf2_crane
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
CRANE_OBJS += crane_model/impl_jac_x.o
CRANE_OBJS += crane_model/impl_jac_xdot.o
CRANE_OBJS += crane_model/impl_jac_u.o

sim_crane_no_interface: $(CRANE_OBJS) sim_crane_no_interface.o
	$(CCC) -o sim_crane_no_interface.out sim_crane_no_interface.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_crane_no_interface build complete."
	@echo

sim_erk_crane: $(CRANE_OBJS) sim_erk_crane.o
	$(CCC) -o sim_erk_crane.out sim_erk_crane.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_erk_crane build complete."
	@echo
	

# WT_OBJS =
# WT_OBJS += simple_wt_model/jac_model.o
# WT_OBJS += simple_wt_model/ode_model.o
# WT_OBJS += simple_wt_model/vde_adj_model.o
# WT_OBJS += simple_wt_model/vde_forw_model.o
# WT_OBJS += simple_wt_model/vde_hess_model.o
# WT_OBJS += sim_erk_wt.o

# sim_erk_wt: $(WT_OBJS)
# 	$(CCC) -o sim_erk_wt.out  $(WT_OBJS) $(LDFLAGS) $(LIBS)
# 	@echo
# 	@echo " Example sim_erk_wt build complete."
# 	@echo


CRANE_GNSF_OBJS =
CRANE_GNSF_OBJS += gnsf_crane_model/res_inc_Jff_fun.o
CRANE_GNSF_OBJS += gnsf_crane_model/get_ints_fun.o
CRANE_GNSF_OBJS += gnsf_crane_model/But_KK_ZZ_LO_fun.o
CRANE_GNSF_OBJS += gnsf_crane_model/f_LO_inc_J_x1k1uz_fun.o
CRANE_GNSF_OBJS += gnsf_crane_model/jac_res_ffx1u_fun.o
CRANE_GNSF_OBJS += sim_gnsf_crane.o

sim_gnsf_crane: $(CRANE_GNSF_OBJS)
	$(CCC) -o sim_gnsf_crane.out  $(CRANE_GNSF_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_gnsf_crane build complete."
	@echo


CRANE_GNSF_EZ_OBJS =
CRANE_GNSF_EZ_OBJS += gnsf_crane_model_empty_z/res_inc_Jff_fun.o
CRANE_GNSF_EZ_OBJS += gnsf_crane_model_empty_z/get_ints_fun.o
CRANE_GNSF_EZ_OBJS += gnsf_crane_model_empty_z/But_KK_ZZ_LO_fun.o
CRANE_GNSF_EZ_OBJS += gnsf_crane_model_empty_z/f_LO_inc_J_x1k1uz_fun.o
CRANE_GNSF_EZ_OBJS += gnsf_crane_model_empty_z/jac_res_ffx1u_fun.o
CRANE_GNSF_EZ_OBJS += sim_gnsf_crane.o # TODO oops

sim_gnsf_crane_empty_z: $(CRANE_GNSF_EZ_OBJS)
	$(CCC) -o sim_gnsf_crane_empty_z.out  $(CRANE_GNSF_EZ_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_gnsf_crane_empty_z build complete."
	@echo

CRANE_GNSF2_OBJS =
CRANE_GNSF2_OBJS += gnsf2_crane_model/But_KK_YY_ZZ_LO_fun.o
CRANE_GNSF2_OBJS += gnsf2_crane_model/f_LO_inc_J_x1k1uz_fun.o
CRANE_GNSF2_OBJS += gnsf2_crane_model/get_ints_fun.o
CRANE_GNSF2_OBJS += gnsf2_crane_model/get_matrices_fun.o
CRANE_GNSF2_OBJS += gnsf2_crane_model/Phi_inc_dy_fun.o
CRANE_GNSF2_OBJS += gnsf2_crane_model/jac_Phi_y_fun.o
CRANE_GNSF2_OBJS += sim_gnsf2_crane.o

sim_gnsf2_crane: $(CRANE_GNSF2_OBJS)
	$(CCC) -o sim_gnsf2_crane.out  $(CRANE_GNSF2_OBJS) $(LDFLAGS) $(LIBS) #/opt/openblas/lib/libopenblas.a -pthread -lgfortran -lm
	@echo
	@echo " Example sim_gnsf2_crane build complete."
	@echo


ANOTHER_CRANE_OBJS =
ANOTHER_CRANE_OBJS += crane_model_2/impl_ode.o
ANOTHER_CRANE_OBJS += crane_model_2/impl_jac_x.o
ANOTHER_CRANE_OBJS += crane_model_2/impl_jac_xdot.o
ANOTHER_CRANE_OBJS += crane_model_2/impl_jac_u.o
sim_irk_crane2: $(ANOTHER_CRANE_OBJS) sim_irk_crane2.o
	$(CCC) -o sim_irk_crane2.out  sim_irk_crane2.o $(ANOTHER_CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_irk_crane2 build complete."
	@echo


BLASFEO_BUG_OBJ =
BLASFEO_BUG_OBJ += blasfeo_bug_minimal_example.o

blasfeo_bug_minimal_example: $(BLASFEO_BUG_OBJ)
	$(CCC) -o blasfeo_bug_minimal_example.out  $(BLASFEO_BUG_OBJ) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example blasfeo_bug_minimal_example build complete."
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


run_sim_crane_example:
	./sim_crane_example.out



#################################################
# wind turbine model with 6 states
#################################################

WT_NX6_OBJS =
WT_NX6_OBJS += wt_model_nx6_expl/expl_ode.o
WT_NX6_OBJS += wt_model_nx6_expl/expl_forw_vde.o

sim_wt_model_nx6: $(WT_NX6_OBJS) sim_wt_model_nx6.o
	$(CCC) -o sim_wt_model_nx6.out sim_wt_model_nx6.o  $(WT_NX6_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_wt_model_nx6 build complete."
	@echo

run_sim_wt_model_nx6:
	./sim_wt_model_nx6.out



WT_NX6P2_OBJS =
WT_NX6P2_OBJS += wt_model_nx6_expl/nx6p2/expl_ode.o
WT_NX6P2_OBJS += wt_model_nx6_expl/nx6p2/expl_forw_vde.o

wind_turbine_nmpc: $(WT_NX6P2_OBJS) wind_turbine_nmpc.o
	$(CCC) -o wind_turbine_nmpc.out wind_turbine_nmpc.o  $(WT_NX6P2_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example wind_turbine_nmpc build complete."
	@echo

run_wind_turbine_nmpc:
	./wind_turbine_nmpc.out



#################################################
# wind turbine model
#################################################

WT_OBJS =
WT_OBJS += wt_model_nx3/ode_energy_balanced_model.o
WT_OBJS += wt_model_nx3/jac_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_adj_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_energy_balanced_model.o
WT_OBJS += wt_model_nx3/vde_hess_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_odeFun_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_jacFun_x_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_jacFun_xdot_energy_balanced_model.o
WT_OBJS += wt_model_nx3/impl_jacFun_u_energy_balanced_model.o

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
CHAIN_OBJS += implicit_chain_model/impl_odeFun_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_x_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_xdot_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_u_chain_nm2.o
# 3 masses
CHAIN_OBJS += implicit_chain_model/impl_odeFun_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_x_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_xdot_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_u_chain_nm3.o
# 4 masses
CHAIN_OBJS += implicit_chain_model/impl_odeFun_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_x_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_xdot_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_u_chain_nm4.o
# 5 masses
CHAIN_OBJS += implicit_chain_model/impl_odeFun_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_x_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_xdot_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_u_chain_nm5.o
# 6 masses
CHAIN_OBJS += implicit_chain_model/impl_odeFun_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_x_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_xdot_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_jacFun_u_chain_nm6.o

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


run_sim_erk_wt:
	./sim_erk_wt.out


#################################################
# clean
#################################################

#################################################
# clean
#################################################

clean:
	rm -f *.o
	rm -f *.out

deep_clean: clean
	rm -f chain_model/*.o
	rm -f wt_model_nx3/*.o
	rm -f wt_model_nx6_expl/*.o
	rm -f wt_model_nx6_expl/nx6p2/*.o

