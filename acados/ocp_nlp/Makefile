# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += ocp_nlp_common.o
OBJS += ocp_nlp_cost_common.o
OBJS += ocp_nlp_cost_ls.o
OBJS += ocp_nlp_cost_nls.o
OBJS += ocp_nlp_cost_external.o
OBJS += ocp_nlp_constraints.o
OBJS += ocp_nlp_dynamics.o
OBJS += ocp_nlp_sqp.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
