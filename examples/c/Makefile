# Author: Gianluca Frison



include ../../Makefile.rule



LIBS =

ifeq ($(ACADOS_WITH_C_INTERFACE), 1)
LIBS += -lacados_c
endif

LIBS = -lacore -lhpipm

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
#EXAMPLES += sim_erk_crane
#EXAMPLES += sim_irk_crane
#EXAMPLES += sim_lifted_irk_crane
#EXAMPLES += mass_spring_example
#EXAMPLES += mass_spring_pcond_hpipm_split
#EXAMPLES += mass_spring_fcond_hpipm_split
#EXAMPLES += mass_spring_offline_fcond_qpoases_split
#EXAMPLES += nonlinear_chain_ocp_nlp
EXAMPLES += sim_crane_no_interface
EXAMPLES += mass_spring_example_no_interface
EXAMPLES += nonlinear_chain_ocp_nlp_no_interface

examples: $(EXAMPLES)



RUN_EXAMPLES =
#RUN_EXAMPLES += run_sim_erk_crane
#RUN_EXAMPLES += run_sim_irk_crane
#RUN_EXAMPLES += run_sim_lifted_irk_crane
#RUN_EXAMPLES += run_mass_spring_example
#RUN_EXAMPLES += run_mass_spring_pcond_hpipm_split
#RUN_EXAMPLES += run_mass_spring_fcond_hpipm_split
#RUN_EXAMPLES += run_mass_spring_offline_fcond_qpoases_split
#RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp
RUN_EXAMPLES += run_sim_crane_no_interface
RUN_EXAMPLES += run_mass_spring_example_no_interface
RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp_no_interface

run_examples: $(RUN_EXAMPLES)



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

sim_irk_crane: $(CRANE_OBJS) sim_irk_crane.o
	$(CCC) -o sim_irk_crane.out  sim_irk_crane.o $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_irk_crane build complete."
	@echo

sim_lifted_irk_crane: $(CRANE_OBJS) sim_lifted_irk_crane.o
	$(CCC) -o sim_lifted_irk_crane.out  sim_lifted_irk_crane.o $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_lifted_irk_crane build complete."
	@echo



MASS_SPRING_OBJS =
MASS_SPRING_OBJS += mass_spring_model/mass_spring_qp.o

mass_spring_example_no_interface: $(MASS_SPRING_OBJS) mass_spring_example_no_interface.o
	$(CCC) -o mass_spring_example_no_interface.out $(MASS_SPRING_OBJS) mass_spring_example_no_interface.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example_no_interface build complete."
	@echo

mass_spring_example: $(MASS_SPRING_OBJS) mass_spring_example.o
	$(CCC) -o mass_spring_example.out $(MASS_SPRING_OBJS) mass_spring_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example build complete."
	@echo

mass_spring_pcond_hpipm_split: mass_spring_pcond_hpipm_split.o
	$(CCC) -o mass_spring_pcond_hpipm_split.out mass_spring_pcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_hpipm_split build complete."
	@echo

mass_spring_fcond_hpipm_split: mass_spring_fcond_hpipm_split.o
	$(CCC) -o mass_spring_fcond_hpipm_split.out mass_spring_fcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_fcond_hpipm_split build complete."
	@echo

mass_spring_offline_fcond_qpoases_split: mass_spring_offline_fcond_qpoases_split.o
	$(CCC) -o mass_spring_offline_fcond_qpoases_split.out mass_spring_offline_fcond_qpoases_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_offline_fcond_qpoases_split build complete."
	@echo



CHAIN_OBJS =
### explicit ###
# 2 masses
CHAIN_OBJS += chain_model/vde_chain_nm2.o
CHAIN_OBJS += chain_model/jac_chain_nm2.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm2.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm2.o
CHAIN_OBJS += chain_model/ls_cost_nm2.o
CHAIN_OBJS += chain_model/ls_costN_nm2.o
# 3 masses
CHAIN_OBJS += chain_model/vde_chain_nm3.o
CHAIN_OBJS += chain_model/jac_chain_nm3.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm3.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm3.o
CHAIN_OBJS += chain_model/ls_cost_nm3.o
CHAIN_OBJS += chain_model/ls_costN_nm3.o
# 4 masses
CHAIN_OBJS += chain_model/vde_chain_nm4.o
CHAIN_OBJS += chain_model/jac_chain_nm4.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm4.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm4.o
CHAIN_OBJS += chain_model/ls_cost_nm4.o
CHAIN_OBJS += chain_model/ls_costN_nm4.o
# 5 masses
CHAIN_OBJS += chain_model/vde_chain_nm5.o
CHAIN_OBJS += chain_model/jac_chain_nm5.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm5.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm5.o
CHAIN_OBJS += chain_model/ls_cost_nm5.o
CHAIN_OBJS += chain_model/ls_costN_nm5.o
# 6 masses
CHAIN_OBJS += chain_model/vde_chain_nm6.o
CHAIN_OBJS += chain_model/jac_chain_nm6.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm6.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm6.o
CHAIN_OBJS += chain_model/ls_cost_nm6.o
CHAIN_OBJS += chain_model/ls_costN_nm6.o
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

nonlinear_chain_ocp_nlp: $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o
	$(CCC) -o nonlinear_chain_ocp_nlp.out $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo



run_sim_crane_no_interface:
	./sim_crane_no_interface.out



run_sim_erk_crane:
	./sim_erk_crane.out



run_sim_irk_crane:
	./sim_irk_crane.out



run_sim_lifted_irk_crane:
	./sim_lifted_irk_crane.out



run_mass_spring_example_no_interface:
	./mass_spring_example_no_interface.out



run_mass_spring_example:
	./mass_spring_example.out



run_mass_spring_pcond_hpipm_split:
	./mass_spring_pcond_hpipm_split.out



run_mass_spring_fcond_hpipm_split:
	./mass_spring_fcond_hpipm_split.out



run_mass_spring_offline_fcond_qpoases_split:
	./mass_spring_offline_fcond_qpoases_split.out



run_nonlinear_chain_ocp_nlp_no_interface:
	./nonlinear_chain_ocp_nlp_no_interface.out



run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out



clean:
	rm -f *.o
	rm -f *.out

deep_clean: clean
	rm -f chain_model/*.o

