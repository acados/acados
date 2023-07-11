#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#




# default make target
all: static_library



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
OBJS += acados/ocp_nlp/ocp_nlp_cost_conl.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_common.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_bgh.o
OBJS += acados/ocp_nlp/ocp_nlp_constraints_bgp.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_common.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_cont.o
OBJS += acados/ocp_nlp/ocp_nlp_dynamics_disc.o
OBJS += acados/ocp_nlp/ocp_nlp_sqp.o
OBJS += acados/ocp_nlp/ocp_nlp_sqp_rti.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_common.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_convexify.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_mirror.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_project.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_project_reduc_hess.o
OBJS += acados/ocp_nlp/ocp_nlp_reg_noreg.o

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
OBJS += acados/ocp_qp/ocp_qp_xcond_solver.o
# sim
OBJS += acados/sim/sim_collocation_utils.o
OBJS += acados/sim/sim_erk_integrator.o
OBJS += acados/sim/sim_irk_integrator.o
OBJS += acados/sim/sim_lifted_irk_integrator.o
OBJS += acados/sim/sim_common.o
OBJS += acados/sim/sim_gnsf.o
# utils
OBJS += acados/utils/math.o
OBJS += acados/utils/print.o
OBJS += acados/utils/timing.o
OBJS += acados/utils/mem.o
OBJS += acados/utils/external_function_generic.o

# C interface
ifeq ($(ACADOS_WITH_C_INTERFACE), 1)
OBJS += interfaces/acados_c/external_function_interface.o
OBJS += interfaces/acados_c/dense_qp_interface.o
OBJS += interfaces/acados_c/ocp_nlp_interface.o
OBJS += interfaces/acados_c/ocp_qp_interface.o
OBJS += interfaces/acados_c/condensing_interface.o
OBJS += interfaces/acados_c/sim_interface.o
endif

# acados dependencies
STATIC_DEPS = blasfeo_static hpipm_static
SHARED_DEPS = blasfeo_shared hpipm_shared
CLEAN_DEPS = blasfeo_clean hpipm_clean
ifeq ($(ACADOS_WITH_QPOASES), 1)
STATIC_DEPS += qpoases_static
SHARED_DEPS += qpoases_shared
CLEAN_DEPS += qpoases_clean
LINK_FLAG_QPOASES = -lqpOASES_e
endif
ifeq ($(ACADOS_WITH_HPMPC), 1)
STATIC_DEPS += hpmpc_static
CLEAN_DEPS += hpmpc_clean
LINK_FLAG_HPMPC = -lhpmpc
endif
ifeq ($(ACADOS_WITH_QPDUNES), 1)
STATIC_DEPS += qpdunes_static
CLEAN_DEPS += qpdunes_clean
LINK_FLAG_QPDUNES = -lqpdunes
endif
ifeq ($(ACADOS_WITH_QORE), 1)
STATIC_DEPS += qore_static
CLEAN_DEPS += qore_clean
LINK_FLAG_QPDUNES = -lqore
endif
ifeq ($(ACADOS_WITH_OSQP), 1)
STATIC_DEPS += osqp_static
SHARED_DEPS += osqp_shared
CLEAN_DEPS += osqp_clean
LINK_FLAG_QPDUNES = -losqp
endif

ifeq ($(ACADOS_WITH_OPENMP), 1)
LINK_FLAG_OPENMP = -fopenmp
endif


static_library: $(STATIC_DEPS)
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	( cd interfaces/acados_c; $(MAKE) obj CC=$(CC) TOP=$(TOP) )
	ar rcs libacados.a $(OBJS)
	mkdir -p lib
	mv libacados.a lib
	mkdir -p include/acados
	cp --parents acados/*/*.h include/
	mkdir -p include/acados_c
	cp -r interfaces/acados_c/*.h include/acados_c
	@echo
	@echo " libacados.a static library build complete."
	@echo

shared_library: link_libs_json $(SHARED_DEPS)
	( cd acados; $(MAKE) obj TOP=$(TOP) )
	( cd interfaces/acados_c; $(MAKE) obj  CC=$(CC) TOP=$(TOP) )
	$(CC) -L./lib -shared -o libacados.so $(OBJS) -lblasfeo -lhpipm -lm -fopenmp
	mkdir -p lib
	mv libacados.so lib
	mkdir -p include/acados
	cp --parents acados/*/*.h include/
	mkdir -p include/acados_c
	cp -r interfaces/acados_c/*.h include/acados_c
	@echo
	@echo " libacados.so shared library build complete."
	@echo

# write linker flags to external libraries into lib/links.json, used in MEX interface
link_libs_json:
	echo "{" > ./lib/link_libs.json
	echo "\t\"openmp\": \"$(LINK_FLAG_OPENMP)\"," >> ./lib/link_libs.json
	echo "\t\"qpoases\": \"$(LINK_FLAG_QPOASES)\"," >> ./lib/link_libs.json
	echo "\t\"qpdunes\": \"$(LINK_FLAG_QPDUNES)\"," >> ./lib/link_libs.json
	echo "\t\"osqp\": \"$(LINK_FLAG_OSQP)\"," >> ./lib/link_libs.json
	echo "\t\"hpmpc\": \"$(LINK_FLAG_HPMPC)\"," >> ./lib/link_libs.json
	echo "\t\"ooqp\": \"$(LINK_FLAG_OOQP)\"" >> ./lib/link_libs.json
	echo "}" >> ./lib/link_libs.json

blasfeo_static:
	( cd $(BLASFEO_PATH); $(MAKE) static_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) MF=PANELMAJ BLAS_API=0 )
	mkdir -p include/blasfeo/include
	mkdir -p lib
	cp $(BLASFEO_PATH)/include/*.h include/blasfeo/include
	cp $(BLASFEO_PATH)/lib/libblasfeo.a lib

blasfeo_shared:
	( cd $(BLASFEO_PATH); $(MAKE) shared_library CC=$(CC) LA=$(BLASFEO_VERSION) TARGET=$(BLASFEO_TARGET) MF=PANELMAJ BLAS_API=0 )
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

qpoases_shared:
	( cd $(QPOASES_PATH); $(MAKE) CC=$(CC) MAKE_STATIC_LIB=0 DLLEXT=so )
	mkdir -p include/qpoases/include
	mkdir -p lib
	cp -r $(QPOASES_PATH)/include/* include/qpoases/include
	cp $(QPOASES_PATH)/bin/libqpOASES_e.so lib

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

osqp_shared: $(OSQP_LIB_SHARED)
	mkdir -p include/osqp/include
	mkdir -p lib
	cp -r $(OSQP_PATH)/include/* include/osqp/include
	mv libosqp.so lib

examples_c: static_library
	( cd examples/c; $(MAKE) examples TOP=$(TOP) )

run_examples_c: examples_c
	( cd examples/c; $(MAKE) run_examples )

unit_tests: static_library
	( cd test; $(MAKE) unit_tests TOP=$(TOP) )

run_unit_tests: unit_tests
	( cd test; $(MAKE) run_unit_tests )

clean:
	( cd acados; $(MAKE) clean )
	( cd examples/c; $(MAKE) clean )
	( cd test; $(MAKE) clean )
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
