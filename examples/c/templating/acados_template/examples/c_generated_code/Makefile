OBJ = 
OBJ+= template_example.o
OBJ+= pendulum_ode_model/pendulum_ode_expl_ode_fun.o
OBJ+= pendulum_ode_model/pendulum_ode_expl_vde_forw.o

INCLUDE_PATH = /usr/local/include
LIB_PATH = ../../../../../../build

all: clean casadi_fun main example

casadi_fun:
	( cd pendulum_ode_model; $(MAKE) pendulum_ode_expl_ode_fun.o pendulum_ode_expl_vde_forw.o )

main:
	gcc -c template_example.c \ 
	-I $(INCLUDE_PATH)/blasfeo/include/ \
	-I $(INCLUDE_PATH)/acados/ -I $(INCLUDE_PATH)/qpOASES_e/ 

example: 
	gcc -o template_example $(OBJ) -L $(LIB_PATH)/acados \
	-L $(LIB_PATH)/external/blasfeo \
	-L $(LIB_PATH)/external/hpipm  \
	-L $(LIB_PATH)/external/qpoases/lib \
	-lacados -lhpipm -lblasfeo -lqpOASES_e -lm \
	-I $(INCLUDE_PATH)/blasfeo/include/ \
	-I $(INCLUDE_PATH)/acados/ \
	-I $(INCLUDE_PATH)/qpOASES_e/  

clean:
	rm -f *.o
	rm -f ./pendulum_ode_model/*.o
