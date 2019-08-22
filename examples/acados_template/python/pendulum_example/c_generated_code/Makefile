#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

ACADOS_FLAGS = -fPIC

ACADOS_FLAGS += -DACADOS_WITH_QPOASES

MODEL_OBJ=

MODEL_OBJ+= pendulum_ode_model/pendulum_ode_expl_ode_fun.o
MODEL_OBJ+= pendulum_ode_model/pendulum_ode_expl_vde_forw.o



OCP_OBJ=


OCP_OBJ+= acados_solver_pendulum_ode.o

SIM_OBJ=
SIM_OBJ+= acados_sim_solver_pendulum_ode.o

EX_OBJ=
EX_OBJ+= main_pendulum_ode.o

OBJ=
OBJ+= $(MODEL_OBJ)
OBJ+= $(SIM_OBJ)
OBJ+= $(OCP_OBJ)



INCLUDE_PATH = /usr/lib/include
LIB_PATH = /usr/lib/lib

all: clean casadi_fun main ocp_solver sim_solver example
shared_lib: bundled_shared_lib ocp_shared_lib sim_shared_lib

CASADI_MODEL_SOURCE=

CASADI_MODEL_SOURCE+= pendulum_ode_expl_ode_fun.c
CASADI_MODEL_SOURCE+= pendulum_ode_expl_vde_forw.c


CASADI_CON_P_SOURCE=

CASADI_CON_H_SOURCE=


casadi_fun:
	( cd pendulum_ode_model; gcc $(ACADOS_FLAGS) -c  $(CASADI_MODEL_SOURCE))
	
	

main:
	gcc $(ACADOS_FLAGS) -c main_pendulum_ode.c -I $(INCLUDE_PATH)/blasfeo/include/ -I $(INCLUDE_PATH)/hpipm/include/ \
	-I $(INCLUDE_PATH) -I $(INCLUDE_PATH)/acados/ -I $(INCLUDE_PATH)/qpOASES_e/

ocp_solver:
	gcc $(ACADOS_FLAGS) -c acados_solver_pendulum_ode.c -I $(INCLUDE_PATH)/blasfeo/include/ -I $(INCLUDE_PATH)/hpipm/include/ \
	-I $(INCLUDE_PATH) -I $(INCLUDE_PATH)/acados/ -I $(INCLUDE_PATH)/qpOASES_e/

sim_solver:
	gcc $(ACADOS_FLAGS) -c acados_sim_solver_pendulum_ode.c -I $(INCLUDE_PATH)/blasfeo/include/ -I $(INCLUDE_PATH)/hpipm/include/ \
	-I $(INCLUDE_PATH) -I $(INCLUDE_PATH)/acados/ -I $(INCLUDE_PATH)/qpOASES_e/

example:
	gcc $(ACADOS_FLAGS) -o main_pendulum_ode $(EX_OBJ) $(OBJ) -L $(LIB_PATH) \
	-L $(LIB_PATH)/acados \
	-L $(LIB_PATH)/external/blasfeo \
	-L $(LIB_PATH)/external/hpipm  \
	-L $(LIB_PATH)/external/qpoases/lib \
	-lacados -lhpipm -lblasfeo -lqpOASES_e -lm \
	-I $(INCLUDE_PATH)/blasfeo/include/ \
	-I $(INCLUDE_PATH)/hpipm/include/ \
	-I $(INCLUDE_PATH) \
	-I $(INCLUDE_PATH)/acados/ \
	-I $(INCLUDE_PATH)/qpOASES_e/

bundled_shared_lib: casadi_fun ocp_solver sim_solver
	gcc $(ACADOS_FLAGS) -shared -o libacados_solver_pendulum_ode.so $(OBJ) -L $(LIB_PATH)/acados \
	-L $(LIB_PATH) \
	-L $(LIB_PATH)/external/blasfeo \
	-L $(LIB_PATH)/external/hpipm  \
	-L $(LIB_PATH)/external/qpoases/lib \
	-lacados -lhpipm -lblasfeo -lqpOASES_e -lm \

ocp_shared_lib: casadi_fun ocp_solver
	gcc $(ACADOS_FLAGS) -shared -o libacados_ocp_solver_pendulum_ode.so $(OCP_OBJ) $(MODEL_OBJ) -L $(LIB_PATH)/acados \
	-L $(LIB_PATH) \
	-L $(LIB_PATH)/external/blasfeo \
	-L $(LIB_PATH)/external/hpipm  \
	-L $(LIB_PATH)/external/qpoases/lib \
	-lacados -lhpipm -lblasfeo -lqpOASES_e -lm \

sim_shared_lib: casadi_fun sim_solver
	gcc $(ACADOS_FLAGS) -shared -o libacados_sim_solver_pendulum_ode.so $(SIM_OBJ) $(MODEL_OBJ) -L $(LIB_PATH)/acados \
	-L $(LIB_PATH) \
	-L $(LIB_PATH)/external/blasfeo \
	-L $(LIB_PATH)/external/hpipm  \
	-L $(LIB_PATH)/external/qpoases/lib \
	-lacados -lhpipm -lblasfeo -lqpOASES_e -lm \

clean:
	rm -f *.o
	rm -f *.so
	rm -f main_pendulum_ode