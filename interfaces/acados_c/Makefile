include ../../Makefile.rule

OBJS =

OBJS += dense_qp/dense_qp_hpipm.o
ifeq ($(ACADOS_WITH_QORE), 1)
OBJS += dense_qp/dense_qp_qore.o
endif
OBJS += dense_qp/dense_qp_qpoases.o
OBJS += ocp_qp/ocp_qp_full_condensing_solver.o
OBJS += ocp_qp/ocp_qp_full_condensing.o
OBJS += ocp_qp/ocp_qp_hpipm.o
ifeq ($(ACADOS_WITH_HPMPC), 1)
OBJS += ocp_qp/ocp_qp_hpmpc.o
endif
OBJS += ocp_qp/ocp_qp_partial_condensing.o
ifeq ($(ACADOS_WITH_QPDUNES), 1)
OBJS += ocp_qp/ocp_qp_qpdunes.o
endif
OBJS += ocp_qp/ocp_qp_sparse_solver.o
OBJS += sim/sim_erk_integrator.o
OBJS += sim/sim_irk_integrator.o
OBJS += sim/sim_lifted_irk_integrator.o
OBJS += utils/casadi_wrapper.o
OBJS += dense_qp.o
OBJS += external_function.o
OBJS += legacy_create.o
OBJS += ocp_lin.o
# OBJS += ocp_nlp.o
OBJS += ocp_qp.o
OBJS += options.o
OBJS += sim.o

static_library: $(OBJS)
	ar rcs libacados_c.a $(OBJS)
	@echo
	@echo " libacados_c.a static library build complete."
	@echo

clean:
	rm -f *.o
	rm -f *.s
	rm -f libacados_c.a
