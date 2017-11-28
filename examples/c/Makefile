# Author: Gianluca Frison

\

include ../../Makefile.rule



examples: mass_spring_hpipm mass_spring_pcond_solver mass_spring_cond_solver mass_spring_pcond_hpipm_split mass_spring_condensing_hpipm_split mass_spring_condensing_hpipm mass_spring_condensing_qpoases mass_spring_condensing_qore nonlinear_chain_ocp_nlp

run_examples: run_mass_spring_hpipm run_mass_spring_pcond_solver run_mass_spring_cond_solver run_mass_spring_pcond_hpipm_split run_mass_spring_condensing_hpipm_split run_mass_spring_condensing_hpipm run_mass_spring_condensing_qpoases run_mass_spring_condensing_qore run_nonlinear_chain_ocp_nlp

mass_spring_hpipm: mass_spring_hpipm.o
	$(CC) -o mass_spring_hpipm.out mass_spring_hpipm.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_hpipm build complete."
	@echo

mass_spring_pcond_solver: mass_spring_pcond_solver.o
	$(CC) -o mass_spring_pcond_solver.out mass_spring_pcond_solver.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_pcond_solver build complete."
	@echo

mass_spring_cond_solver: mass_spring_cond_solver.o
	$(CC) -o mass_spring_cond_solver.out mass_spring_cond_solver.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_cond_solver build complete."
	@echo

mass_spring_pcond_hpipm_split: mass_spring_pcond_hpipm_split.o
	$(CC) -o mass_spring_pcond_hpipm_split.out mass_spring_pcond_hpipm_split.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_pcond_hpipm_split build complete."
	@echo

mass_spring_condensing_hpipm_split: mass_spring_condensing_hpipm_split.o
	$(CC) -o mass_spring_condensing_hpipm_split.out mass_spring_condensing_hpipm_split.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_condensing_hpipm_split build complete."
	@echo

mass_spring_condensing_hpipm: mass_spring_condensing_hpipm.o
	$(CC) -o mass_spring_condensing_hpipm.out mass_spring_condensing_hpipm.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_condensing_hpipm build complete."
	@echo

mass_spring_condensing_qpoases: mass_spring_condensing_qpoases.o
	$(CC) -o mass_spring_condensing_qpoases.out mass_spring_condensing_qpoases.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example mass_spring_condensing_qpoases build complete."
	@echo

mass_spring_condensing_qore: mass_spring_condensing_qore.o
	$(CC) -o mass_spring_condensing_qore.out mass_spring_condensing_qore.o $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
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
	$(CC) -o nonlinear_chain_ocp_nlp.out $(OBJS) $(LDFLAGS) -lacore -lhpipm -lqore_dense -lblasfeo -lqpOASES_e -lm -lblas
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo

run_mass_spring_hpipm:
	./mass_spring_hpipm.out

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
