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



# define sources and use make's implicit rules to generate object files (*.o)

# model
MODEL_SRC=
MODEL_SRC+= pendulum_ode_model/pendulum_ode_impl_dae_fun.c
MODEL_SRC+= pendulum_ode_model/pendulum_ode_impl_dae_fun_jac_x_xdot_z.c
MODEL_SRC+= pendulum_ode_model/pendulum_ode_impl_dae_jac_x_xdot_u_z.c
MODEL_OBJ := $(MODEL_SRC:.c=.o)
# optimal control problem - mostly CasADi exports
OCP_SRC=
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_0_fun.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_0_fun_jac_ut_xt.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_0_hess.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_fun.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_fun_jac_ut_xt.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_hess.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_e_fun.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_e_fun_jac_ut_xt.c
OCP_SRC+= pendulum_ode_cost/pendulum_ode_cost_y_e_hess.c

OCP_SRC+= acados_solver_pendulum_ode.c
OCP_OBJ := $(OCP_SRC:.c=.o)

# for target example
EX_SRC= main_pendulum_ode.c
EX_OBJ := $(EX_SRC:.c=.o)
EX_EXE := $(EX_SRC:.c=)


# combine model, (potentially) sim and ocp object files
OBJ=
OBJ+= $(MODEL_OBJ)
OBJ+= $(OCP_OBJ)

EXTERNAL_DIR=
EXTERNAL_LIB=

INCLUDE_PATH = ${ACADOS_INSTALL_DIR}/include
LIB_PATH = ${ACADOS_INSTALL_DIR}/lib

# preprocessor flags for make's implicit rules
CPPFLAGS+= -I$(INCLUDE_PATH)
CPPFLAGS+= -I$(INCLUDE_PATH)/acados
CPPFLAGS+= -I$(INCLUDE_PATH)/blasfeo/include
CPPFLAGS+= -I$(INCLUDE_PATH)/hpipm/include


# define the c-compiler flags for make's implicit rules
CFLAGS = -fPIC -std=c99   -O2#-fno-diagnostics-show-line-numbers -g

# # Debugging
# CFLAGS += -g3

# linker flags
LDFLAGS+= -L$(LIB_PATH)


# link to libraries
LDLIBS+= -lacados
LDLIBS+= -lhpipm
LDLIBS+= -lblasfeo
LDLIBS+= -lm

# libraries
LIBACADOS_SOLVER=libacados_solver_pendulum_ode.so
LIBACADOS_OCP_SOLVER=libacados_ocp_solver_pendulum_ode.so
LIBACADOS_SIM_SOLVER=lib$(SIM_SRC:.c=.so)

# virtual targets
.PHONY : all clean

all: clean example
shared_lib: bundled_shared_lib ocp_shared_lib sim_shared_lib

# some linker targets
example: $(EX_OBJ) $(OBJ)
	$(CC) $^ -o $(EX_EXE) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(OBJ) $(EX_OBJ) $(EX_SIM_OBJ)
	$(RM) $(LIBACADOS_SOLVER) $(LIBACADOS_OCP_SOLVER) $(LIBACADOS_SIM_SOLVER)
	$(RM) $(EX_EXE) $(EX_SIM_EXE)
