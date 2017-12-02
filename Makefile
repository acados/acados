# Author: Gianluca Frison



include ./Makefile.rule



OBJS =

# ocp nlp
OBJS += acados/ocp_nlp/ocp_nlp_common.o
# OBJS += acados/ocp_nlp/ocp_nlp_sm_gn.o
OBJS += acados/ocp_nlp/ocp_nlp_gn_sqp.o
# dense qp
OBJS += acados/dense_qp/dense_qp_common.o
OBJS += acados/dense_qp/dense_qp_hpipm.o
OBJS += acados/dense_qp/dense_qp_qpoases.o
OBJS += acados/dense_qp/dense_qp_qore.o
# ocp qp
OBJS += acados/ocp_qp/ocp_qp_common.o
OBJS += acados/ocp_qp/ocp_qp_common_frontend.o
OBJS += acados/ocp_qp/ocp_qp_hpipm.o
OBJS += acados/ocp_qp/ocp_qp_hpmpc.o
OBJS += acados/ocp_qp/ocp_qp_condensing.o
OBJS += acados/ocp_qp/ocp_qp_partial_condensing.o
OBJS += acados/ocp_qp/ocp_qp_sparse_solver.o
OBJS += acados/ocp_qp/ocp_qp_condensing_solver.o
#sim
OBJS += acados/sim/sim_casadi_wrapper.o
OBJS += acados/sim/sim_collocation.o
OBJS += acados/sim/sim_erk_integrator.o
OBJS += acados/sim/sim_lifted_irk_integrator.o
OBJS += acados/sim/sim_common.o
# utils
OBJS += acados/utils/math.o
OBJS += acados/utils/copy.o
OBJS += acados/utils/print.o
OBJS += acados/utils/timing.o
OBJS += acados/utils/mem.o
OBJS += acados/utils/create.o


static_library: blasfeo_static hpipm_static hpmpc_static qpoases_static qore_static
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

hpmpc_static: blasfeo_static
	( cd external/hpmpc; $(MAKE) static_library CC=$(CC) TARGET=$(HPMPC_TARGET) BLASFEO_PATH=$(TOP)/external/blasfeo  )
	mkdir -p include/hpmpc
	mkdir -p lib
	cp external/hpmpc/include/*.h include/hpmpc
	cp external/hpmpc/libhpmpc.a lib

qpoases_static:
	( cd external/qpoases; $(MAKE) CC=$(CC) )
	mkdir -p include/qpoases
	mkdir -p lib
	cp -r external/qpoases/include/* include/qpoases
	cp external/qpoases/bin/libqpOASES_e.a lib

qore_static: blasfeo_static
	mkdir -p external/qore/external/blasfeo
	cp external/blasfeo/include/*.h external/qore/external/blasfeo
	cp external/blasfeo/lib/libblasfeo.a external/qore/external/blasfeo
	( cd external/qore; $(MAKE) static_dense; )
	mkdir -p include/qore
	mkdir -p lib
	cp external/qore/qp_types.h include/qore
	#cp external/qore/KKTPACK_DENSE/include/*.h include/qore
	#cp external/qore/KKTPACK_DENSE/source/*.h include/qore
	#cp external/qore/QPCORE/include/*.h include/qore
	cp external/qore/QPSOLVER_DENSE/include/*.h include/qore
	#cp external/qore/QPSOLVER_DENSE/source/*.h include/qore
	cp external/qore/bin/libqore_dense.a lib

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
	( cd external/hpmpc; $(MAKE) clean )
	( cd external/qpoases; $(MAKE) clean )
	( cd external/qore; $(MAKE) purge )
	rm -rf include
	rm -rf lib
