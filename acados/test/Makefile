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

