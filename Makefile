# Author: Gianluca Frison



include ./Makefile.rule



OBJS =

# ocp nlp
OBJS += acados/ocp_nlp/allocate_ocp_nlp.o
OBJS += acados/ocp_nlp/ocp_nlp_common.o
# OBJS += acados/ocp_nlp/ocp_nlp_common_ext_dep.o
# OBJS += acados/ocp_nlp/ocp_nlp_sm_gn.o
OBJS += acados/ocp_nlp/ocp_nlp_gn_sqp.o
# OBJS += acados/ocp_nlp/ocp_nlp_sqp_ext_dep.o
# dense qp
OBJS += acados/dense_qp/dense_qp_common.o
OBJS += acados/dense_qp/dense_qp_common_ext_dep.o
OBJS += acados/dense_qp/dense_qp_hpipm.o
OBJS += acados/dense_qp/dense_qp_hpipm_ext_dep.o
OBJS += acados/dense_qp/dense_qp_qpoases.o
OBJS += acados/dense_qp/dense_qp_qpoases_ext_dep.o
# ocp qp
OBJS += acados/ocp_qp/ocp_qp_common.o
OBJS += acados/ocp_qp/ocp_qp_common_ext_dep.o
OBJS += acados/ocp_qp/ocp_qp_common_frontend.o
OBJS += acados/ocp_qp/ocp_qp_hpipm.o
OBJS += acados/ocp_qp/ocp_qp_hpipm_ext_dep.o
OBJS += acados/ocp_qp/ocp_qp_condensing_hpipm.o
OBJS += acados/ocp_qp/ocp_qp_condensing_hpipm_ext_dep.o
OBJS += acados/ocp_qp/ocp_qp_condensing_qpoases.o
OBJS += acados/ocp_qp/ocp_qp_condensing_qpoases_ext_dep.o
OBJS += acados/ocp_qp/ocp_qp_condensing.o
OBJS += acados/ocp_qp/ocp_qp_condensing_ext_dep.o
#sim
OBJS += acados/sim/allocate_sim.o
OBJS += acados/sim/casadi_wrapper.o
OBJS += acados/sim/sim_collocation.o
OBJS += acados/sim/sim_erk_integrator.o
OBJS += acados/sim/sim_lifted_irk_integrator.o
# utils
OBJS += acados/utils/math.o
OBJS += acados/utils/copy.o
OBJS += acados/utils/print.o
OBJS += acados/utils/timing.o
OBJS += acados/utils/mem.o


static_library:
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	ar rcs libacore.a $(OBJS)
	mkdir -p lib
	mv libacore.a lib
	@echo
	@echo " libacore.a static library build complete."
	@echo

examples_c:
	( cd examples/c; $(MAKE) examples TOP=$(TOP) )

run_examples_c:
	( cd examples/c; $(MAKE) run_examples )

clean:
	( cd acados; $(MAKE) clean )
	( cd examples/c; $(MAKE) clean )
	rm -f libacore.a
	rm -rf lib
