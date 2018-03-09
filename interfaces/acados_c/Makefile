include ../../Makefile.rule

OBJS =

OBJS += dense_qp/dense_qp_hpipm_i.o
ifeq ($(ACADOS_WITH_QORE), 1)
OBJS += dense_qp/dense_qp_qore_i.o
endif
OBJS += dense_qp/dense_qp_qpoases_i.o
OBJS += ocp_qp/ocp_qp_full_condensing_solver_i.o
OBJS += ocp_qp/ocp_qp_full_condensing_i.o
OBJS += ocp_qp/ocp_qp_hpipm_i.o
ifeq ($(ACADOS_WITH_HPMPC), 1)
OBJS += ocp_qp/ocp_qp_hpmpc_i.o
endif
OBJS += ocp_qp/ocp_qp_partial_condensing_i.o
ifeq ($(ACADOS_WITH_QPDUNES), 1)
OBJS += ocp_qp/ocp_qp_qpdunes_i.o
endif
OBJS += ocp_qp/ocp_qp_partial_condensing_solver_i.o
OBJS += sim/sim_erk_integrator_i.o
OBJS += sim/sim_lifted_irk_integrator_i.o
OBJS += dense_qp_i.o
OBJS += legacy_create_i.o
# OBJS += ocp_lin_i.o
# OBJS += ocp_nlp.o
OBJS += ocp_qp_interface.o
OBJS += options_i.o
OBJS += sim_i.o
OBJS += external_function_generic.o

static_library: $(OBJS)
	ar rcs libacados_c.a $(OBJS)
	@echo
	@echo " libacados_c.a static library build complete."
	@echo

clean:
	rm -f *.o
	rm -f *.s
	rm -f libacados_c.a
