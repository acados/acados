# Author: Gianluca Frison



# default make target
all: acados_static acados_c_static



# include config & tailored rules
include ./Makefile.rule
include ./Makefile.osqp



# acados sources
OBJS =

# ocp nlp
OBJS += acados/ocp_nlp/ocp_nlp_common.o
OBJS += acados/ocp_nlp/ocp_nlp_cost_common.o
OBJS += acados/ocp_nlp/ocp_nlp_cost_ls.o
OBJS += acados/ocp_nlp/ocp_nlp_cost_nls.o
OBJS += acados/ocp_nlp/ocp_nlp_cost_external.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_common.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_bgh.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_bghp.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_common.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_cont.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_disc.o
OBJS += acados/ocp_nlp/ocp_nlp_sqp.o
OBJS += acados/ocp_nlp/ocp_nlp_sqp_rti.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_common.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_conv.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_mirror.o

# dense qp
OBJS += acados/dense_qp/dense_qp_common.o
OBJS += acados/dense_qp/dense_qp_hpipm.o
ifeq ($(ACADOS_WITH_QPOASES), 1)
OBJS += acados/dense_qp/dense_qp_qpoases.o
endif
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
ifeq ($(ACADOS_WITH_OSQP), 1)
OBJS += acados/ocp_qp/ocp_qp_osqp.o
endif
OBJS += acados/ocp_qp/ocp_qp_partial_condensing.o
OBJS += acados/ocp_qp/ocp_qp_full_condensing.o
OBJS += acados/ocp_qp/ocp_qp_partial_condensing_solver.o
OBJS += acados/ocp_qp/ocp_qp_full_condensing_solver.o
# sim
OBJS += acados/sim/sim_collocation_utils.o
OBJS += acados/sim/sim_erk_integrator.o
OBJS += acados/sim/sim_irk_integrator.o
OBJS += acados/sim/sim_lifted_irk_integrator.o
OBJS += acados/sim/sim_common.o
OBJS += acados/sim/sim_gnsf.o
# utils
OBJS += acados/utils/math.o
OBJS += acados/utils/copy.o
OBJS += acados/utils/print.o
OBJS += acados/utils/timing.o
OBJS += acados/utils/mem.o
OBJS += acados/utils/external_function_generic.o



# acados dependencies
STATIC_DEPS = blasfeo_static hpipm_static
SHARED_DEPS = blasfeo_shared hpipm_shared
CLEAN_DEPS = blasfeo_clean hpipm_clean
ifeq ($(ACADOS_WITH_QPOASES), 1)
STATIC_DEPS += qpoases_static
CLEAN_DEPS += qpoases_clean
endif
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
ifeq ($(ACADOS_WITH_OSQP), 1)
STATIC_DEPS += osqp_static
CLEAN_DEPS += osqp_clean
endif



acados_static: $(STATIC_DEPS)
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	ar rcs libacore.a $(OBJS)
	mkdir -p lib
	mv libacore.a lib
	@echo
	@echo " libacore.a static library build complete."
	@echo

acados_shared: $(SHARED_DEPS)
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	ar rcs libacore.a $(OBJS)
	$(CC) -L./lib -shared -o libacore.so $(OBJS) -lblasfeo -lhpipm -lm -fopenmp
	mkdir -p lib
	mv libacore.so lib
	@echo
	@echo " libacore.so shared library build complete."
	@echo

blasfeo_static:
	( cd $(BLASFEO_PATH); $(MAKE) static_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) BLAS_API=0 )
	mkdir -p include/blasfeo/include
	mkdir -p lib
	cp $(BLASFEO_PATH)/include/*.h include/blasfeo/include
	cp $(BLASFEO_PATH)/lib/libblasfeo.a lib

blasfeo_shared:
	( cd $(BLASFEO_PATH); $(MAKE) shared_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) BLAS_API=0 )
	mkdir -p include/blasfeo/include
	mkdir -p lib
	cp $(BLASFEO_PATH)/include/*.h include/blasfeo/include
	cp $(BLASFEO_PATH)/lib/libblasfeo.so lib

hpipm_static: blasfeo_static
	( cd $(HPIPM_PATH); $(MAKE) static_library CC=$(CC) TARGET=$(HPIPM_TARGET) BLASFEO_PATH=$(BLASFEO_PATH) )
	mkdir -p include/hpipm/include
	mkdir -p lib
	cp $(HPIPM_PATH)/include/*.h include/hpipm/include
	cp $(HPIPM_PATH)/lib/libhpipm.a lib

hpipm_shared: blasfeo_shared
	( cd $(HPIPM_PATH); $(MAKE) shared_library CC=$(CC) TARGET=$(HPIPM_TARGET) BLASFEO_PATH=$(BLASFEO_PATH) )
	mkdir -p include/hpipm/include
	mkdir -p lib
	cp $(HPIPM_PATH)/include/*.h include/hpipm/include
	cp $(HPIPM_PATH)/lib/libhpipm.so lib

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
	mkdir -p include/qore/QPSOLVER_DENSE/include
	mkdir -p include/qore/QPSOLVER_DENSE/source
	mkdir -p include/qore/KKTPACK_DENSE/include
	mkdir -p include/qore/KKTPACK_DENSE/source
	mkdir -p include/qore/QORE/include
	mkdir -p lib
	cp $(QORE_PATH)/qp_types.h include/qore/
	cp $(QORE_PATH)/QPSOLVER_DENSE/include/*.h include/qore/QPSOLVER_DENSE/include
	cp $(QORE_PATH)/QPSOLVER_DENSE/source/*.h include/qore/QPSOLVER_DENSE/source
	cp $(QORE_PATH)/KKTPACK_DENSE/source/*.h include/qore/KKTPACK_DENSE/source
	cp $(QORE_PATH)/KKTPACK_DENSE/include/*.h include/qore/KKTPACK_DENSE/include
	cp $(QORE_PATH)/QPCORE/include/*.h include/qore/QORE/include
	cp $(QORE_PATH)/bin/libqore_dense.a lib

qpdunes_static:
	( cd $(QPDUNES_PATH); $(MAKE) CC=$(CC) )
	mkdir -p include/qpdunes/include
	mkdir -p lib
	cp -r $(QPDUNES_PATH)/include/* include/qpdunes/include
	cp $(QPDUNES_PATH)/src/libqpdunes.a lib
	cp $(QPDUNES_PATH)/externals/qpOASES-3.0beta/bin/libqpOASES.a lib
	
osqp_static: $(OSQP_LIB_STATIC)
	mkdir -p include/osqp/include
	mkdir -p lib
	cp -r $(OSQP_PATH)/include/* include/osqp/include
	mv libosqp.a lib

acados_c_static: acados_static
ifeq ($(ACADOS_WITH_C_INTERFACE), 1)
	( cd interfaces/acados_c; $(MAKE) static_library CC=$(CC) TOP=$(TOP) )
	mkdir -p include/acados_c
	mkdir -p lib
	cp -r interfaces/acados_c/*.h include/acados_c
	mv interfaces/acados_c/libacados_c.a lib
endif

acados_c_shared: acados_shared
ifeq ($(ACADOS_WITH_C_INTERFACE), 1)
	( cd interfaces/acados_c; $(MAKE) shared_library CC=$(CC) TOP=$(TOP) )
	mkdir -p include/acados_c
	mkdir -p lib
	cp -r interfaces/acados_c/*.h include/acados_c
	mv interfaces/acados_c/libacados_c.so lib
endif

examples_c: acados_c_static
	( cd examples/c; $(MAKE) examples TOP=$(TOP) )

run_examples_c: examples_c
	( cd examples/c; $(MAKE) run_examples )

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
	
osqp_clean:
	@$(RM) $(OSQP_ALL_OBJ)
	@$(RM) $(OSQP_QDLDL_INC_DIR)qdldl_types.h
	@$(RM) $(OSQP_INC_DIR)osqp_configure.h

deep_clean: clean $(CLEAN_DEPS)
	( cd examples/c; $(MAKE) clean )

clean_models:
	( cd examples/c; $(MAKE) clean_models )

purge: deep_clean clean_models
	rm -rf include
	rm -rf lib
