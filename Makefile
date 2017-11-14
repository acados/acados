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
OBJS += acados/ocp_qp/ocp_qp_partial_condensing.o
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
OBJS += acados/utils/create.o


static_library: blasfeo_static hpipm_static qpoases_static
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	ar rcs libacore.a $(OBJS)
	mkdir -p lib
	mv libacore.a lib
	@echo
	@echo " libacore.a static library build complete."
	@echo

blasfeo_static:
	( cd external/blasfeo; $(MAKE) static_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) )
	mkdir -p include/blasfeo
	mkdir -p lib
	cp external/blasfeo/include/*.h include/blasfeo
	cp external/blasfeo/lib/libblasfeo.a lib

hpipm_static: blasfeo_static
	( cd external/hpipm; $(MAKE) static_library CC=$(CC) TARGET=$(HPIPM_TARGET) BLASFEO_PATH=$(TOP)/external/blasfeo )
	mkdir -p include/hpipm
	mkdir -p lib
	cp external/hpipm/include/*.h include/hpipm
	cp external/hpipm/lib/libhpipm.a lib

qpoases_static:
	( cd external/qpoases; $(MAKE) CC=$(CC) )
	mkdir -p include/qpoases
	mkdir -p lib
	cp -r external/qpoases/include/* include/qpoases
	cp external/qpoases/bin/libqpOASES_e.a lib

examples_c:
	( cd examples/c; $(MAKE) examples TOP=$(TOP) )

run_examples_c:
	( cd examples/c; $(MAKE) run_examples )

run_example_chain:
	( cd examples/c; $(MAKE) run_nonlinear_chain_ocp_nlp )

clean:
	( cd acados; $(MAKE) clean )
	( cd examples/c; $(MAKE) clean )

deep_clean: clean
	( cd external/blasfeo; $(MAKE) clean )
	( cd external/hpipm; $(MAKE) clean )
	( cd external/qpoases; $(MAKE) clean )
	rm -rf include
	rm -rf lib
