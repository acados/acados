# Author: Gianluca Frison



include ./Makefile.rule



# acados sources
OBJS =

# ocp nlp
OBJS += acados/ocp_nlp/ocp_nlp_common.o
# OBJS += acados/ocp_nlp/ocp_nlp_sm_gn.o
# OBJS += acados/ocp_nlp/ocp_nlp_gn_sqp.o
# dense qp
OBJS += acados/dense_qp/dense_qp_common.o
OBJS += acados/dense_qp/dense_qp_hpipm.o
OBJS += acados/dense_qp/dense_qp_qpoases.o
ifeq ($(ACADOS_WITH_QORE), 1)
OBJS += acados/dense_qp/dense_qp_qore.o
endif
# ocp qp
OBJS += acados/ocp_qp/ocp_qp_common.o
OBJS += acados/ocp_qp/ocp_qp_common_frontend.o
OBJS += acados/ocp_qp/ocp_qp_hpipm.o
ifeq ($(ACADOS_WITH_HPMPC), 1)
OBJS += acados/ocp_qp/ocp_qp_hpmpc.o
endif
ifeq ($(ACADOS_WITH_QPDUNES), 1)
OBJS += acados/ocp_qp/ocp_qp_qpdunes.o
endif
OBJS += acados/ocp_qp/ocp_qp_partial_condensing.o
OBJS += acados/ocp_qp/ocp_qp_full_condensing.o
OBJS += acados/ocp_qp/ocp_qp_sparse_solver.o
OBJS += acados/ocp_qp/ocp_qp_full_condensing_solver.o
# sim
OBJS += acados/sim/sim_collocation.o
OBJS += acados/sim/sim_erk_integrator.o
OBJS += acados/sim/sim_lifted_irk_integrator.o
OBJS += acados/sim/sim_common.o
# utils
OBJS += acados/utils/casadi_wrapper.o
OBJS += acados/utils/math.o
OBJS += acados/utils/copy.o
OBJS += acados/utils/external_function.o
OBJS += acados/utils/print.o
OBJS += acados/utils/timing.o
OBJS += acados/utils/mem.o
OBJS += acados/utils/external_function_generic.o



# acados dependencies
STATIC_DEPS = blasfeo_static hpipm_static qpoases_static
CLEAN_DEPS = blasfeo_clean hpipm_clean qpoases_clean
ifeq ($(ACADOS_WITH_HPMPC), 1)
STATIC_DEPS += hpmpc_static
CLEAN_DEPS += hpmpc_clean
endif
ifeq ($(ACADOS_WITH_QPDUNES), 1)
STATIC_DEPS += qpdunes_static
CLEAN_DEPS += qpdunes_clean
endif
ifeq ($(ACADOS_WITH_QORE), 1)
STATIC_DEPS += qore_static
CLEAN_DEPS += qore_clean
endif



all: acados_c_static

acados_static: $(STATIC_DEPS)
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	ar rcs libacore.a $(OBJS)
	mkdir -p lib
	mv libacore.a lib
	@echo
	@echo " libacore.a static library build complete."
	@echo

blasfeo_static:
	( cd $(BLASFEO_PATH); $(MAKE) static_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) )
	mkdir -p include/blasfeo/include
	mkdir -p lib
	cp $(BLASFEO_PATH)/include/*.h include/blasfeo/include
	cp $(BLASFEO_PATH)/lib/libblasfeo.a lib

hpipm_static: blasfeo_static
	( cd $(HPIPM_PATH); $(MAKE) static_library CC=$(CC) TARGET=$(HPIPM_TARGET) BLASFEO_PATH=$(BLASFEO_PATH) )
	mkdir -p include/hpipm/include
	mkdir -p lib
	cp $(HPIPM_PATH)/include/*.h include/hpipm/include
	cp $(HPIPM_PATH)/lib/libhpipm.a lib

hpmpc_static: blasfeo_static
	( cd $(HPMPC_PATH); $(MAKE) static_library CC=$(CC) TARGET=$(HPMPC_TARGET) BLASFEO_PATH=$(BLASFEO_PATH)  )
	mkdir -p include/hpmpc/include
	mkdir -p lib
	cp $(HPMPC_PATH)/include/*.h include/hpmpc/include
	cp $(HPMPC_PATH)/libhpmpc.a lib

qpoases_static:
	( cd $(QPOASES_PATH); $(MAKE) CC=$(CC) )
	mkdir -p include/qpoases/include
	mkdir -p lib
	cp -r $(QPOASES_PATH)/include/* include/qpoases/include
	cp $(QPOASES_PATH)/bin/libqpOASES_e.a lib

# TODO how is BLASFEO path set for QORE ?????
qore_static: blasfeo_static
	( cd $(QORE_PATH); $(MAKE) static_dense; )
	mkdir -p include/qore/include
	mkdir -p lib
	cp $(QORE_PATH)/qp_types.h include/qore/include
	cp $(QORE_PATH)/QPSOLVER_DENSE/include/*.h include/qore/include
	cp $(QORE_PATH)/QPSOLVER_DENSE/source/*.h include/qore/include
	cp $(QORE_PATH)/KKTPACK_DENSE/source/*.h include/qore/include
	cp $(QORE_PATH)/KKTPACK_DENSE/include/*.h include/qore/include
	cp $(QORE_PATH)/QPCORE/include/*.h include/qore/include
	cp $(QORE_PATH)/bin/libqore_dense.a lib

qpdunes_static:
	( cd $(QPDUNES_PATH); $(MAKE) CC=$(CC) )
	mkdir -p include/qpdunes/include
	mkdir -p lib
	cp -r $(QPDUNES_PATH)/include/* include/qpdunes/include
	cp $(QPDUNES_PATH)/src/libqpdunes.a lib
	cp $(QPDUNES_PATH)/externals/qpOASES-3.0beta/bin/libqpOASES.a lib

acados_c_static: acados_static
	( cd interfaces/acados_c; $(MAKE) static_library CC=$(CC) TOP=$(TOP) )
	mkdir -p include/acados_c
	mkdir -p include/acados_c/dense_qp
	# mkdir -p include/acados_c/ocp_lin
	# mkdir -p include/acados_c/ocp_nlp
	mkdir -p include/acados_c/ocp_qp
	mkdir -p include/acados_c/sim
	mkdir -p lib
	cp -r interfaces/acados_c/*.h include/acados_c
	cp -r interfaces/acados_c/dense_qp/*.h include/acados_c/dense_qp
	# cp -r interfaces/acados_c/ocp_lin/*.h include/acados_c/ocp_lin
	# cp -r interfaces/acados_c/ocp_nlp/*.h include/acados_c/ocp_nlp
	cp -r interfaces/acados_c/ocp_qp/*.h include/acados_c/ocp_qp
	cp -r interfaces/acados_c/sim/*.h include/acados_c/sim
	mv interfaces/acados_c/libacados_c.a lib

examples_c: acados_c_static
	( cd examples/c; $(MAKE) examples TOP=$(TOP) )

run_examples_c: examples_c
	( cd examples/c; $(MAKE) run_examples )

run_example_chain:
	( cd examples/c; $(MAKE) run_nonlinear_chain_ocp_nlp )

clean:
	( cd acados; $(MAKE) clean )
	( cd examples/c; $(MAKE) clean )
	( cd interfaces/acados_c; $(MAKE) clean )

blasfeo_clean:
	( cd $(BLASFEO_PATH); $(MAKE) clean )

hpipm_clean:
	( cd $(HPIPM_PATH); $(MAKE) clean )

hpmpc_clean:
	( cd $(HPMPC_PATH); $(MAKE) clean )

qpoases_clean:
	( cd $(QPOASES_PATH); $(MAKE) clean )

qore_clean:
	( cd $(QORE_PATH); $(MAKE) purge )

qpdunes_clean:
	( cd $(QPDUNES_PATH); $(MAKE) clean )

deep_clean: clean $(CLEAN_DEPS)
	( cd examples/c; $(MAKE) deep_clean )
	rm -rf include
	rm -rf lib
