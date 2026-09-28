#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.




include ../Makefile.rule

CPPFLAGS += -I$(TOP)/external



LIBS =

LIBS += -lacados -lhpipm

ifeq ($(ACADOS_WITH_QPOASES), 1)
LIBS += -lqpOASES_e
endif

ifeq ($(ACADOS_WITH_HPMPC), 1)
LIBS += -lhpmpc
endif

ifeq ($(ACADOS_WITH_QPDUNES), 1)
CCC = $(CXX)
LIBS += -lqpdunes -lqpOASES
else
CCC = $(CC)
endif

ifeq ($(ACADOS_WITH_QORE), 1)
LIBS += -lqore_dense
endif

ifeq ($(ACADOS_WITH_OSQP), 1)
LIBS += -losqp -ldl
endif

LIBS += -lblasfeo -lm -lblas -llapack

ifeq ($(ACADOS_WITH_OPENMP), 1)
LIBS += -fopenmp
endif



TESTS =
TESTS += test_qpsolvers

unit_tests: $(TESTS)



RUN_TESTS =
RUN_TESTS += run_test_qpsolvers

run_unit_tests: $(RUN_UNIT_TESTS)



#################################################
# (linear) mass spring model
#################################################

MASS_SPRING_OBJS =
MASS_SPRING_OBJS += ../examples/c/no_interface_examples/mass_spring_model/mass_spring_qp.o



test_qpsolvers: $(MASS_SPRING_OBJS) ocp_qp/test_qpsolvers.o
	$(CCC) -o test_qpsolvers.out $(MASS_SPRING_OBJS) ocp_qp/test_qpsolvers.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example test_qpsolvers build complete."
	@echo

run_test_qpsolvers:
	./test_qpsolvers.out



#################################################
# clean
#################################################

clean:
	rm -f *.o
	rm -f *.out

