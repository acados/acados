#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


include ../../Makefile.rule

OBJS =

OBJS += ocp_qp_common.o
OBJS += ocp_qp_common_frontend.o
OBJS += ocp_qp_hpipm.o
ifeq ($(ACADOS_WITH_HPMPC), 1)
OBJS += ocp_qp_hpmpc.o
endif
ifeq ($(ACADOS_WITH_QPDUNES), 1)
OBJS += ocp_qp_qpdunes.o
endif
ifeq ($(ACADOS_WITH_OSQP), 1)
OBJS += ocp_qp_osqp.o
endif
OBJS += ocp_qp_partial_condensing.o
OBJS += ocp_qp_full_condensing.o
OBJS += ocp_qp_xcond_solver.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
