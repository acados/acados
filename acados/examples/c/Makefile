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




include ../../Makefile.rule



LIBS =

LIBS += -lacados -lhpipm

ifeq ($(ACADOS_WITH_QPOASES), 1)
LIBS += -lqpOASES_e
endif

ifeq ($(ACADOS_WITH_HPMPC), 1)
LIBS += -lhpmpc
endif

ifeq ($(ACADOS_WITH_QPDUNES), 1)
CCC = $(CXX)
LIBS += -lqpdunes -lqpOASES
else
CCC = $(CC)
endif

ifeq ($(ACADOS_WITH_QORE), 1)
LIBS += -lqore_dense
endif

ifeq ($(ACADOS_WITH_OSQP), 1)
LIBS += -losqp -ldl
endif

LIBS += -lblasfeo -lm -lblas -llapack

ifeq ($(ACADOS_WITH_OPENMP), 1)
LIBS += -fopenmp
endif


# Comment this out to enable using gprof
# CFLAGS  += -pg
# LDFLAGS += -pg

EXAMPLES =
EXAMPLES += dense_qp
EXAMPLES += ocp_qp
EXAMPLES += sim_wt_model_nx3
EXAMPLES += sim_wt_model_nx6
EXAMPLES += sim_pendulum_dae
EXAMPLES += sim_crane_example
EXAMPLES += sim_gnsf_crane
EXAMPLES += mass_spring_example
EXAMPLES += mass_spring_nmpc_example
##EXAMPLES += mass_spring_pcond_split
##EXAMPLES += mass_spring_fcond_split
##EXAMPLES += mass_spring_offline_fcond_qpoases_split
EXAMPLES += nonlinear_chain_ocp_nlp
# EXAMPLES += sim_crane_no_interface
##EXAMPLES += mass_spring_example_no_interface
#EXAMPLES += nonlinear_chain_ocp_nlp_no_interface
EXAMPLES += wind_turbine_nmpc
#EXAMPLES += engine_example
EXAMPLES += regularization
EXAMPLES += simple_dae_example

examples: $(EXAMPLES)


# run all
RUN_EXAMPLES =
RUN_EXAMPLES += run_dense_qp
RUN_EXAMPLES += run_ocp_qp
RUN_EXAMPLES += run_sim_wt_model_nx3
RUN_EXAMPLES += run_sim_wt_model_nx6
RUN_EXAMPLES += run_sim_pendulum_dae
RUN_EXAMPLES += run_sim_crane_example
RUN_EXAMPLES += run_sim_gnsf_crane
RUN_EXAMPLES += run_mass_spring_example
RUN_EXAMPLES += run_mass_spring_nmpc_example
##RUN_EXAMPLES += run_mass_spring_pcond_split
##RUN_EXAMPLES += run_mass_spring_fcond_split
#RUN_EXAMPLES += run_mass_spring_offline_fcond_qpoases_split
RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp
# RUN_EXAMPLES += run_sim_crane_no_interface
##RUN_EXAMPLES += run_mass_spring_example_no_interface
#RUN_EXAMPLES += run_nonlinear_chain_ocp_nlp_no_interface
RUN_EXAMPLES += run_wind_turbine_nmpc
#RUN_EXAMPLES += run_engine_example
RUN_EXAMPLES += run_regularization
RUN_EXAMPLES += run_simple_dae_example

# overwrite to run just one
#RUN_EXAMPLES = run_wind_turbine_nmpc
#RUN_EXAMPLES = run_mass_spring_nmpc_example
#RUN_EXAMPLES = run_regularization
#RUN_EXAMPLES = run_mass_spring_example
#RUN_EXAMPLES = run_nonlinear_chain_ocp_nlp

run_examples: $(RUN_EXAMPLES)

#################################################
# crane model
#################################################

CRANE_OBJS =
CRANE_OBJS += crane_model/ode_model.o
CRANE_OBJS += crane_model/jac_model.o
CRANE_OBJS += crane_model/vde_adj_model.o
CRANE_OBJS += crane_model/vde_forw_model.o
CRANE_OBJS += crane_model/vde_hess_model.o
CRANE_OBJS += crane_model/impl_ode_fun.o
CRANE_OBJS += crane_model/impl_ode_fun_jac_x_xdot.o
CRANE_OBJS += crane_model/impl_ode_jac_x_xdot_u.o

#sim_crane_no_interface: $(CRANE_OBJS) sim_crane_no_interface.o
#	$(CCC) -o sim_crane_no_interface.out sim_crane_no_interface.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
#	@echo
#	@echo " Example sim_crane_no_interface build complete."
#	@echo

#run_sim_crane_no_interface:
#	./sim_crane_no_interface.out

#sim_erk_crane: $(CRANE_OBJS) sim_erk_crane.o
#	$(CCC) -o sim_erk_crane.out sim_erk_crane.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
#	@echo
#	@echo " Example sim_erk_crane build complete."
#	@echo

sim_crane_example: $(CRANE_OBJS) sim_crane_example.o
	$(CCC) -o sim_crane_example.out sim_crane_example.o  $(CRANE_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_crane_example build complete."
	@echo

run_sim_crane_example:
	./sim_crane_example.out


CRANE_GNSF_OBJS =
CRANE_GNSF_OBJS += crane_nx9_model/crane_nx9_phi_fun.o
CRANE_GNSF_OBJS += crane_nx9_model/crane_nx9_phi_fun_jac_y.o
CRANE_GNSF_OBJS += crane_nx9_model/crane_nx9_f_lo_fun_jac_x1k1uz.o
CRANE_GNSF_OBJS += crane_nx9_model/crane_nx9_get_matrices_fun.o
CRANE_GNSF_OBJS += crane_nx9_model/crane_nx9_phi_jac_y_uhat.o
CRANE_GNSF_OBJS += sim_gnsf_crane.o

sim_gnsf_crane: $(CRANE_GNSF_OBJS)
	$(CCC) -o sim_gnsf_crane.out  $(CRANE_GNSF_OBJS) $(LDFLAGS) $(LIBS) #/opt/openblas/lib/libopenblas.a -pthread -lgfortran -lm
	@echo
	@echo " Example sim_gnsf_crane build complete."
	@echo

run_sim_gnsf_crane:
	./sim_gnsf_crane.out



#################################################
# wind turbine model with 6 states
#################################################

WT_NX6_OBJS =
WT_NX6_OBJS += wt_model_nx6/expl_ode_fun.o
WT_NX6_OBJS += wt_model_nx6/expl_vde_for.o
WT_NX6_OBJS += wt_model_nx6/expl_vde_adj.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_fun.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_fun_jac_x_xdot.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_fun_jac_x_xdot_u.o
WT_NX6_OBJS += wt_model_nx6/impl_ode_jac_x_xdot_u.o
WT_NX6_OBJS += wt_model_nx6/phi_fun.o
WT_NX6_OBJS += wt_model_nx6/phi_fun_jac_y.o
WT_NX6_OBJS += wt_model_nx6/phi_jac_y_uhat.o
WT_NX6_OBJS += wt_model_nx6/f_lo_fun_jac_x1k1uz.o
WT_NX6_OBJS += wt_model_nx6/get_matrices_fun.o

sim_wt_model_nx6: $(WT_NX6_OBJS) sim_wt_model_nx6.o
	$(CCC) -o sim_wt_model_nx6.out sim_wt_model_nx6.o  $(WT_NX6_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_wt_model_nx6 build complete."
	@echo

run_sim_wt_model_nx6:
	./sim_wt_model_nx6.out



WT_NX6P2_OBJS =
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_expl_ode_fun.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_expl_vde_for.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_fun_jac_x_xdot_u.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_impl_ode_jac_x_xdot_u.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_phi_fun.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_phi_fun_jac_y.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_phi_jac_y_uhat.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_f_lo_fun_jac_x1k1uz.o
WT_NX6P2_OBJS += wt_model_nx6/nx6p2/wt_nx6p2_get_matrices_fun.o

wind_turbine_nmpc: $(WT_NX6P2_OBJS) wind_turbine_nmpc.o
	$(CCC) -o wind_turbine_nmpc.out wind_turbine_nmpc.o  $(WT_NX6P2_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example wind_turbine_nmpc build complete."
	@echo

run_wind_turbine_nmpc:
	./wind_turbine_nmpc.out


#################################################
# wind turbine model with 3 states
#################################################

WT_OBJS =
WT_OBJS += wt_model_nx3/expl_ode_fun.o
WT_OBJS += wt_model_nx3/expl_vde_for.o
WT_OBJS += wt_model_nx3/expl_vde_adj.o
WT_OBJS += wt_model_nx3/impl_ode_fun.o
WT_OBJS += wt_model_nx3/impl_ode_fun_jac_x_xdot.o
WT_OBJS += wt_model_nx3/impl_ode_jac_x_xdot_u.o
WT_OBJS += wt_model_nx3/impl_ode_fun_jac_x_xdot_u.o
WT_OBJS += wt_model_nx3/phi_fun.o
WT_OBJS += wt_model_nx3/phi_fun_jac_y.o
WT_OBJS += wt_model_nx3/phi_jac_y_uhat.o
WT_OBJS += wt_model_nx3/f_lo_fun_jac_x1k1uz.o
WT_OBJS += wt_model_nx3/get_matrices_fun.o

sim_wt_model_nx3: $(WT_OBJS) sim_wt_model_nx3.o
	$(CCC) -o sim_wt_model_nx3.out sim_wt_model_nx3.o  $(WT_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_wt_model_nx3 build complete."
	@echo

run_sim_wt_model_nx3:
	./sim_wt_model_nx3.out

#################################################
# inverted_pendulum model (DAE)
#################################################

INV_PENDULUM_OBJS =
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun_jac_x_xdot.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_impl_ode_jac_x_xdot_u.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_impl_ode_fun_jac_x_xdot_u.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_impl_ode_hess.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_fun.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_fun_jac_y.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_gnsf_phi_jac_y_uhat.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_gnsf_f_lo_fun_jac_x1k1uz.o
INV_PENDULUM_OBJS += pendulum_dae_model/pendulum_dae_dyn_gnsf_get_matrices_fun.o

sim_pendulum_dae: $(INV_PENDULUM_OBJS) sim_pendulum_dae.o
	$(CCC) -o sim_pendulum_dae.out sim_pendulum_dae.o  $(INV_PENDULUM_OBJS) $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example sim_pendulum_dae build complete."
	@echo

run_sim_pendulum_dae:
	./sim_pendulum_dae.out



#################################################
# (linear) mass spring model
#################################################

MASS_SPRING_OBJS =
MASS_SPRING_OBJS += no_interface_examples/mass_spring_model/mass_spring_qp.o

mass_spring_example_no_interface: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_example_no_interface.o
	$(CCC) -o mass_spring_example_no_interface.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_example_no_interface.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example_no_interface build complete."
	@echo

run_mass_spring_example_no_interface:
	./mass_spring_example_no_interface.out



mass_spring_example: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_example.o
	$(CCC) -o mass_spring_example.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_example build complete."
	@echo

run_mass_spring_example:
	./mass_spring_example.out



mass_spring_nmpc_example: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_nmpc_example.o
	$(CCC) -o mass_spring_nmpc_example.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_nmpc_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_nmpc_example build complete."
	@echo

run_mass_spring_nmpc_example:
	./mass_spring_nmpc_example.out



mass_spring_pcond_split: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_pcond_split.o
	$(CCC) -o mass_spring_pcond_split.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_pcond_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_pcond_split build complete."
	@echo

run_mass_spring_pcond_split:
	./mass_spring_pcond_split.out



mass_spring_fcond_split: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_fcond_split.o
	$(CCC) -o mass_spring_fcond_split.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_fcond_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_fcond_split build complete."
	@echo

run_mass_spring_fcond_split:
	./mass_spring_fcond_split.out



mass_spring_offline_fcond_qpoases_split: $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_offline_fcond_qpoases_split.o
	$(CCC) -o mass_spring_offline_fcond_qpoases_split.out $(MASS_SPRING_OBJS) no_interface_examples/mass_spring_offline_fcond_qpoases_split.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example mass_spring_offline_fcond_qpoases_split build complete."
	@echo

run_mass_spring_offline_fcond_qpoases_split:
	./mass_spring_offline_fcond_qpoases_split.out



#################################################
# (nonlinear) chain model
#################################################

CHAIN_OBJS =
### explicit ###
# 2 masses
CHAIN_OBJS += chain_model/vde_chain_nm2.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm2.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm2.o
CHAIN_OBJS += chain_model/ls_cost_nm2.o
CHAIN_OBJS += chain_model/ls_costN_nm2.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm2.o
# 3 masses
CHAIN_OBJS += chain_model/vde_chain_nm3.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm3.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm3.o
CHAIN_OBJS += chain_model/ls_cost_nm3.o
CHAIN_OBJS += chain_model/ls_costN_nm3.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm3.o
# 4 masses
CHAIN_OBJS += chain_model/vde_chain_nm4.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm4.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm4.o
CHAIN_OBJS += chain_model/ls_cost_nm4.o
CHAIN_OBJS += chain_model/ls_costN_nm4.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm4.o
# 5 masses
CHAIN_OBJS += chain_model/vde_chain_nm5.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm5.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm5.o
CHAIN_OBJS += chain_model/ls_cost_nm5.o
CHAIN_OBJS += chain_model/ls_costN_nm5.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm5.o
# 6 masses
CHAIN_OBJS += chain_model/vde_chain_nm6.o
CHAIN_OBJS += chain_model/vde_adj_chain_nm6.o
#CHAIN_OBJS += chain_model/vde_hess_chain_nm6.o
CHAIN_OBJS += chain_model/ls_cost_nm6.o
CHAIN_OBJS += chain_model/ls_costN_nm6.o
CHAIN_OBJS += chain_model/casadi_erk4_chain_nm6.o
# external cost
CHAIN_OBJS += chain_model/chain_nm_2_external_cost.o
CHAIN_OBJS += chain_model/chain_nm_3_external_cost.o
CHAIN_OBJS += chain_model/chain_nm_4_external_cost.o
CHAIN_OBJS += chain_model/chain_nm_5_external_cost.o
CHAIN_OBJS += chain_model/chain_nm_6_external_cost.o
### implicit ###
# 2 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm2.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm2.o
# 3 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm3.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm3.o
# 4 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm4.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm4.o
# 5 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm5.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm5.o
# 6 masses
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_fun_jac_x_xdot_u_chain_nm6.o
CHAIN_OBJS += implicit_chain_model/impl_ode_jac_x_xdot_u_chain_nm6.o


nonlinear_chain_ocp_nlp_no_interface: $(CHAIN_OBJS) no_interface_examples/nonlinear_chain_ocp_nlp_no_interface.o
	$(CCC) -o nonlinear_chain_ocp_nlp_no_interface.out $(CHAIN_OBJS) no_interface_examples/nonlinear_chain_ocp_nlp_no_interface.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp_no_interface build complete."
	@echo

run_nonlinear_chain_ocp_nlp_no_interface:
	./nonlinear_chain_ocp_nlp_no_interface.out



nonlinear_chain_ocp_nlp: $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o
	$(CCC) -o nonlinear_chain_ocp_nlp.out $(CHAIN_OBJS) nonlinear_chain_ocp_nlp.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example nonlinear_chain_ocp_nlp build complete."
	@echo

run_nonlinear_chain_ocp_nlp:
	./nonlinear_chain_ocp_nlp.out



#################################################
# dense qp
#################################################

DENSE_QP = 

dense_qp: $(DENSE_QP) dense_qp.o
	$(CXX) -o dense_qp.out $(DENSE_QP) dense_qp.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example dense_qp build complete."
	@echo

run_dense_qp:
	./dense_qp.out




#################################################
# ocp qp
#################################################

OCP_QP = 

ocp_qp: $(OCP_QP) ocp_qp.o
	$(CCC) -o ocp_qp.out $(OCP_QP) ocp_qp.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example ocp_qp build complete."
	@echo

run_ocp_qp:
	./ocp_qp.out




#################################################
# pendulum model
#################################################

ENGINE_OBJS = 
ENGINE_OBJS += engine_model/engine_impl_dae_fun.o
ENGINE_OBJS += engine_model/engine_impl_dae_fun_jac_x_xdot_z.o
ENGINE_OBJS += engine_model/engine_impl_dae_jac_x_xdot_u_z.o
ENGINE_OBJS += engine_model/engine_ls_cost.o
ENGINE_OBJS += engine_model/engine_ls_cost_N.o

engine_example: $(ENGINE_OBJS) engine_example.o
	$(CXX) -o engine_example.out $(ENGINE_OBJS) engine_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example engine_example build complete."
	@echo

run_engine_example:
	./engine_example.out




#################################################
# regularization
#################################################

REGULARIZATION = 

regularization: $(REGULARIZATION) regularization.o
	$(CCC) -o regularization.out $(REGULARIZATION) regularization.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example regularization build complete."
	@echo

run_regularization:
	./regularization.out




#################################################
# simple dae example
#################################################

SIMPLE_DAE_OBJS =
SIMPLE_DAE_OBJS += simple_dae_model/simple_dae_impl_ode_fun.o
SIMPLE_DAE_OBJS += simple_dae_model/simple_dae_impl_ode_jac_x_xdot_u_z.o
SIMPLE_DAE_OBJS += simple_dae_model/simple_dae_impl_ode_fun_jac_x_xdot_z.o

simple_dae_example: $(SIMPLE_DAE_OBJS) simple_dae_example.o
	$(CCC) -o simple_dae_example.out $(SIMPLE_DAE_OBJS) simple_dae_example.o $(LDFLAGS) $(LIBS)
	@echo
	@echo " Example simple_dae_example build complete."
	@echo

run_simple_dae_example:
	./simple_dae_example.out




#################################################
# clean
#################################################

clean:
	rm -f *.o
	rm -f *.out
	rm -f no_interface_examples/*.o
	rm -f no_interface_examples/*.out

clean_models:
	rm -f chain_model/*.o
	rm -f implicit_chain_model/*.o
	rm -f crane_model/*.o
	rm -f crane_nx9_model/*.o
	rm -f wt_model_nx3/*.o
	rm -f wt_model_nx6/*.o
	rm -f wt_model_nx6/nx6p2/*.o
	rm -f pendulum_dae_model/*.o
	rm -f pendulum_model/*.o

