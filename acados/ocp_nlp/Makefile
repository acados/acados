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


include ../../Makefile.rule

OBJS =

OBJS += ocp_nlp_common.o
OBJS += ocp_nlp_cost_common.o
OBJS += ocp_nlp_cost_ls.o
OBJS += ocp_nlp_cost_nls.o
OBJS += ocp_nlp_cost_external.o
OBJS += ocp_nlp_constraints_common.o
OBJS += ocp_nlp_constraints_bgh.o
OBJS += ocp_nlp_constraints_bgp.o
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
