# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += ocp_nlp_common.o
OBJS += ocp_nlp_cost_common.o
OBJS += ocp_nlp_cost_ls.o
OBJS += ocp_nlp_cost_nls.o
OBJS += ocp_nlp_cost_external.o
OBJS += ocp_nlp_constraints_common.o
OBJS += ocp_nlp_constraints_bgh.o
OBJS += ocp_nlp_constraints_bghp.o
OBJS += ocp_nlp_dynamics_common.o
OBJS += ocp_nlp_dynamics_cont.o
OBJS += ocp_nlp_dynamics_disc.o
OBJS += ocp_nlp_sqp.o
OBJS += ocp_nlp_sqp_rti.o
OBJS += ocp_nlp_reg_common.o
OBJS += ocp_nlp_reg_convexify.o
OBJS += ocp_nlp_reg_mirror.o
OBJS += ocp_nlp_reg_project.o
OBJS += ocp_nlp_reg_project_reduc_hess.o
OBJS += ocp_nlp_reg_noreg.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
