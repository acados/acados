# Author: Gianluca Frison

include ../../Makefile.rule

ifeq ($(USE_QPDUNES), 1)
CCC = $(CXX)
LIBS = -lacados_c -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e  -lqpdunes -lqpOASES -lm -lblas
else
CCC = $(CC)
LIBS = -lacados_c -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
endif

examples: mass_spring_example sim_erk_crane mass_spring_pcond_hpipm_split mass_spring_fcond_hpipm_split nonlinear_chain_ocp_nlp

run_examples: run_mass_spring_example run_sim_erk_crane run_mass_spring_pcond_hpipm_split run_mass_spring_fcond_hpipm_split run_nonlinear_chain_ocp_nlp

CRANE_OBJS =
CRANE_OBJS += crane_model/jac_model.o
CRANE_OBJS += crane_model/ode_model.o
CRANE_OBJS += crane_model/vde_adj_model.o
CRANE_OBJS += crane_model/vde_forw_model.o
CRANE_OBJS += crane_model/vde_hess_model.o
CRANE_OBJS += sim_erk_crane.o

sim_erk_crane: $(CRANE_OBJS)
	$(CCC) -o sim_erk_crane.out  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_erk_crane build complete."
	@echo

# sim_irk_crane_yt: $(CRANE_OBJS)
# 	$(CCC) -o sim_irk_crane_yt.out  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
# 	@echo
# 	@echo " Example sim_irk_crane_yt build complete."
# 	@echo

mass_spring_example: mass_spring_example.o
	$(CCC) -o mass_spring_example.out mass_spring_example.o $(LDFLAGS) $(LIBS)
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

OBJS =
OBJS += chain_model/vde_chain_nm2.o
OBJS += chain_model/jac_chain_nm2.o
OBJS += chain_model/vde_adj_chain_nm2.o
OBJS += chain_model/vde_hess_chain_nm2.o
OBJS += chain_model/vde_chain_nm3.o
OBJS += chain_model/jac_chain_nm3.o
OBJS += chain_model/vde_adj_chain_nm3.o
OBJS += chain_model/vde_hess_chain_nm3.o
OBJS += chain_model/vde_chain_nm4.o
OBJS += chain_model/jac_chain_nm4.o
OBJS += chain_model/vde_adj_chain_nm4.o
OBJS += chain_model/vde_hess_chain_nm4.o
OBJS += nonlinear_chain_ocp_nlp.o

nonlinear_chain_ocp_nlp: $(OBJS)
	$(CCC) -o nonlinear_chain_ocp_nlp.out $(OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo

run_sim_erk_crane:
	./sim_erk_crane.out

run_mass_spring_example:
	./mass_spring_example.out

run_mass_spring_pcond_hpipm_split:
	./mass_spring_pcond_hpipm_split.out

run_mass_spring_fcond_hpipm_split:
	./mass_spring_fcond_hpipm_split.out

run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out

clean:
	rm -f chain_model/*.o
	rm -f *.o
	rm -f *.out
