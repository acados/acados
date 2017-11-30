# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += ocp_qp_common.o
OBJS += ocp_qp_common_frontend.o
OBJS += ocp_qp_hpipm.o
OBJS += ocp_qp_hpmpc.o
OBJS += ocp_qp_condensing.o
OBJS += ocp_qp_partial_condensing.o
OBJS += ocp_qp_sparse_solver.o
OBJS += ocp_qp_condensing_solver.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
