# Author Gianluca Frison

include ../../Makefile.rule

OBJS =

OBJS += allocate_sim.o
OBJS += sim_casadi_wrapper.o
OBJS += sim_collocation.o
OBJS += sim_erk_integrator.o
OBJS += sim_lifted_irk_integrator.o
OBJS += sim_erk_integrator_yt.o
OBJS += sim_rk_common_yt.o
OBJS += sim_common_yt.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s
