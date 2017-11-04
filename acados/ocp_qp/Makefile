# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += ocp_qp_common.o
# OBJS += ocp_qp_common_ext_dep.o
# OBJS += ocp_qp_common_frontend.o
OBJS += ocp_qp_hpipm.o
# OBJS += ocp_qp_hpipm_ext_dep.o
OBJS += ocp_qp_condensing_hpipm.o
# OBJS += ocp_qp_condensing_hpipm_ext_dep.o
OBJS += ocp_qp_condensing_qpoases.o
# OBJS += ocp_qp_condensing_qpoases_ext_dep.o
# OBJS += ocp_qp_condensing.o
# OBJS += ocp_qp_condensing_ext_dep.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
