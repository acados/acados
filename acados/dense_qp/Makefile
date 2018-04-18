# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += dense_qp_common.o
OBJS += dense_qp_hpipm.o
ifeq ($(ACADOS_WITH_QPOASES), 1)
OBJS += dense_qp_qpoases.o
endif
ifeq ($(ACADOS_WITH_QORE), 1)
OBJS += dense_qp_qore.o
endif

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
