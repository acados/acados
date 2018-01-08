# Author: Gianluca Frison

\

include ../../Makefile.rule


ifeq ($(USE_QPDUNES), 1)
CCC = $(CXX)
LIBS = -lacore -lhpipm -lhpmpc -lqore_dense -lblasfeo -lqpOASES_e  -lqpdunes -lqpOASES -lm -lblas
else
CCC = $(CC)
LIBS = -lacore -lhpipm -lhpmpc -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
endif

examples: sim_erk_crane mass_spring_hpipm mass_spring_hpmpc mass_spring_qpdunes mass_spring_pcond_solver mass_spring_cond_solver mass_spring_pcond_hpipm_split mass_spring_condensing_hpipm_split mass_spring_condensing_hpipm mass_spring_condensing_qpoases mass_spring_condensing_qore nonlinear_chain_ocp_nlp

run_examples: run_sim_erk_crane run_mass_spring_hpipm run_mass_spring_hpmpc run_mass_spring_qpdunes run_mass_spring_pcond_solver run_mass_spring_cond_solver run_mass_spring_pcond_hpipm_split run_mass_spring_condensing_hpipm_split run_mass_spring_condensing_hpipm run_mass_spring_condensing_qpoases run_mass_spring_condensing_qore run_nonlinear_chain_ocp_nlp

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

sim_irk_crane_yt: $(CRANE_OBJS)
	$(CCC) -o sim_irk_crane_yt.out  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_irk_crane_yt build complete."
	@echo

mass_spring_hpipm: mass_spring_hpipm.o
	$(CCC) -o mass_spring_hpipm.out mass_spring_hpipm.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_hpipm build complete."
	@echo

mass_spring_hpmpc: mass_spring_hpmpc_lazarus.o
	$(CCC) -o mass_spring_hpmpc_lazarus.out mass_spring_hpmpc_lazarus.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_hpmpc build complete."
	@echo

mass_spring_qpdunes: mass_spring_qpdunes.o
	$(CCC) -o mass_spring_qpdunes.out mass_spring_qpdunes.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_qpdunes build complete."
	@echo

mass_spring_pcond_solver: mass_spring_pcond_solver.o
	$(CCC) -o mass_spring_pcond_solver.out mass_spring_pcond_solver.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_solver build complete."
	@echo

mass_spring_cond_solver: mass_spring_cond_solver.o
	$(CCC) -o mass_spring_cond_solver.out mass_spring_cond_solver.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_cond_solver build complete."
	@echo

mass_spring_pcond_hpipm_split: mass_spring_pcond_hpipm_split.o
	$(CCC) -o mass_spring_pcond_hpipm_split.out mass_spring_pcond_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_hpipm_split build complete."
	@echo

mass_spring_condensing_hpipm_split: mass_spring_condensing_hpipm_split.o
	$(CCC) -o mass_spring_condensing_hpipm_split.out mass_spring_condensing_hpipm_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_condensing_hpipm_split build complete."
	@echo

mass_spring_condensing_hpipm: mass_spring_condensing_hpipm.o
	$(CCC) -o mass_spring_condensing_hpipm.out mass_spring_condensing_hpipm.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_condensing_hpipm build complete."
	@echo

mass_spring_condensing_qpoases: mass_spring_condensing_qpoases.o
	$(CCC) -o mass_spring_condensing_qpoases.out mass_spring_condensing_qpoases.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_condensing_qpoases build complete."
	@echo

mass_spring_condensing_qore: mass_spring_condensing_qore.o
	$(CCC) -o mass_spring_condensing_qore.out mass_spring_condensing_qore.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_condensing_qore build complete."
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

run_mass_spring_hpipm:
	./mass_spring_hpipm.out

run_mass_spring_hpmpc:
	./mass_spring_hpmpc_lazarus.out

run_mass_spring_qpdunes:
	./mass_spring_qpdunes.out

run_mass_spring_pcond_hpipm_split:
	./mass_spring_pcond_hpipm_split.out

run_mass_spring_pcond_solver:
	./mass_spring_pcond_solver.out

run_mass_spring_cond_solver:
	./mass_spring_cond_solver.out

run_mass_spring_condensing_hpipm_split:
	./mass_spring_condensing_hpipm_split.out

run_mass_spring_condensing_hpipm:
	./mass_spring_condensing_hpipm.out

run_mass_spring_condensing_qpoases:
	./mass_spring_condensing_qpoases.out

run_mass_spring_condensing_qore:
	./mass_spring_condensing_qore.out

run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out

clean:
	rm -f chain_model/*.o
	rm -f *.o
	rm -f *.out
