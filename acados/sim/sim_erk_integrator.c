/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/mem.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"

#include "acados/utils/external_function.h"



#define FW_VDE_NUMIN 5
#define FW_VDE_NUMOUT 3
#define ADJ_VDE_NUMIN 4
#define ADJ_VDE_NUMOUT 1
#define HESS_VDE_NUMIN 6
#define HESS_VDE_NUMOUT 2



int sim_erk_forward_vde_input_dims[FW_VDE_NUMIN] = {0};
int sim_erk_forward_vde_output_dims[FW_VDE_NUMOUT] = {0};
external_function_dims sim_erk_forward_vde_dims = {
    FW_VDE_NUMIN,
    FW_VDE_NUMOUT,
    sim_erk_forward_vde_input_dims,
    sim_erk_forward_vde_output_dims
};



int sim_erk_adjoint_vde_input_dims[ADJ_VDE_NUMIN] = {0};
int sim_erk_adjoint_vde_output_dims[ADJ_VDE_NUMOUT] = {0};
external_function_dims sim_erk_adjoint_vde_dims = {
    ADJ_VDE_NUMIN,
    ADJ_VDE_NUMOUT,
    sim_erk_adjoint_vde_input_dims,
    sim_erk_adjoint_vde_output_dims
};



int sim_erk_hess_vde_input_dims[HESS_VDE_NUMIN] = {0};
int sim_erk_hess_vde_output_dims[HESS_VDE_NUMOUT] = {0};
external_function_dims sim_erk_hess_vde_dims = {
    HESS_VDE_NUMIN,
    HESS_VDE_NUMOUT,
    sim_erk_hess_vde_input_dims,
    sim_erk_hess_vde_output_dims
};



int sim_erk_integrator_calculate_args_size(sim_dims *dims, void *submodules_)
{
    sim_erk_integrator_submodules *submodules = (sim_erk_integrator_submodules *) submodules_;
 
    int size = sizeof(sim_erk_integrator_args);

    int ns = dims->num_stages;
    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    size += 3*sizeof(external_function_fcn_ptrs);

    if (submodules->forward_vde != NULL) {
        size += submodules->forward_vde->calculate_args_size(&sim_erk_forward_vde_dims, submodules->forward_vde->submodules);
    }
    
    if (submodules->adjoint_vde != NULL) {
        size += submodules->adjoint_vde->calculate_args_size(&sim_erk_adjoint_vde_dims, submodules->adjoint_vde->submodules);
    }
    
    if (submodules->hess_vde != NULL) {
        size += submodules->hess_vde->calculate_args_size(&sim_erk_hess_vde_dims, submodules->hess_vde->submodules);
    }
    
    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_erk_integrator_assign_args(sim_dims *dims, void **submodules_, void *raw_memory)
{
    sim_erk_integrator_submodules *submodules = (sim_erk_integrator_submodules *) *submodules_;
    
    char *c_ptr = (char *) raw_memory;

    sim_erk_integrator_args *args = (sim_erk_integrator_args *) c_ptr;
    c_ptr += sizeof(sim_erk_integrator_args);

    int ns = dims->num_stages;
    args->num_stages = ns;

    align_char_to(8, &c_ptr);

    assign_double(ns*ns, &args->A_mat, &c_ptr);
    assign_double(ns, &args->b_vec, &c_ptr);
    assign_double(ns, &args->c_vec, &c_ptr);

    if (submodules->forward_vde != NULL) {
        args->submodules.forward_vde = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *forward_vde_submodules = submodules->forward_vde->submodules;
        args->forward_vde_args = submodules->forward_vde->assign_args(&sim_erk_forward_vde_dims, &(forward_vde_submodules), c_ptr);
        c_ptr += submodules->forward_vde->calculate_args_size(&sim_erk_forward_vde_dims, submodules->forward_vde->submodules);

        *(args->submodules.forward_vde) = *(submodules->forward_vde);
        args->submodules.forward_vde->submodules = forward_vde_submodules;
    } else {
        args->submodules.forward_vde = NULL;
    }

    if (submodules->adjoint_vde != NULL) {
        args->submodules.adjoint_vde = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *adjoint_vde_submodules = submodules->adjoint_vde->submodules;
        args->adjoint_vde_args = submodules->adjoint_vde->assign_args(&sim_erk_adjoint_vde_dims, &(adjoint_vde_submodules), c_ptr);
        c_ptr += submodules->adjoint_vde->calculate_args_size(&sim_erk_adjoint_vde_dims, submodules->adjoint_vde->submodules);

        *(args->submodules.adjoint_vde) = *(submodules->adjoint_vde);
        args->submodules.adjoint_vde->submodules = adjoint_vde_submodules;
    } else {
        args->submodules.adjoint_vde = NULL;
    }

    if (submodules->hess_vde != NULL) {
        args->submodules.hess_vde = (external_function_fcn_ptrs *)c_ptr;
        c_ptr += sizeof(external_function_fcn_ptrs);

        void *hess_vde_submodules = submodules->hess_vde->submodules;
        args->hess_vde_args = submodules->hess_vde->assign_args(&sim_erk_hess_vde_dims, &(hess_vde_submodules), c_ptr);
        c_ptr += submodules->hess_vde->calculate_args_size(&sim_erk_hess_vde_dims, submodules->hess_vde->submodules);

        *(args->submodules.hess_vde) = *(submodules->hess_vde);
        args->submodules.hess_vde->submodules = hess_vde_submodules;
    } else {
        args->submodules.hess_vde = NULL;
    }

    assert((char*)raw_memory + sim_erk_integrator_calculate_args_size(dims, *submodules_) >= c_ptr);

    // Update submodules pointer
    *submodules_ = (void *) &(args->submodules);

    return (void *)args;
}



void *sim_erk_integrator_copy_args(sim_dims *dims, void *raw_memory, void *source_)
{
    sim_erk_integrator_args *source = (sim_erk_integrator_args *)source_;
    sim_erk_integrator_args *dest;

    sim_erk_integrator_submodules *submodules = &source->submodules;

    dest = sim_erk_integrator_assign_args(dims, (void **) &submodules, raw_memory);

    dest->interval = source->interval;
    dest->num_stages = source->num_stages;
    dest->num_steps = source->num_steps;
    dest->num_forw_sens = source->num_forw_sens;
    dest->sens_forw = source->sens_forw;
    dest->sens_adj = source->sens_adj;
    dest->sens_hess = source->sens_hess;

    int ns = dims->num_stages;

    memcpy(dest->A_mat, source->A_mat, ns*ns*sizeof(double));
    memcpy(dest->c_vec, source->c_vec, ns*sizeof(double));
    memcpy(dest->b_vec, source->b_vec, ns*sizeof(double));

    source->submodules.forward_vde->copy_args(&sim_erk_forward_vde_dims, dest->forward_vde_args, source->forward_vde_args);
    source->submodules.adjoint_vde->copy_args(&sim_erk_adjoint_vde_dims, dest->adjoint_vde_args, source->adjoint_vde_args);
    source->submodules.hess_vde->copy_args(&sim_erk_hess_vde_dims, dest->hess_vde_args, source->hess_vde_args);

    return (void *)dest;
}



void sim_erk_integrator_initialize_default_args(sim_dims *dims, void *args_)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *) args_;
    int ns = args->num_stages;

    assert(args->num_stages == 4 && "only number of stages = 4 implemented!");

    memcpy(args->A_mat,((double[]){0, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1, 0, 0, 0, 0}),
        sizeof(*args->A_mat) * (ns * ns));
    memcpy(args->b_vec, ((double[]){1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6}),
        sizeof(*args->b_vec) * (ns));
    memcpy(args->c_vec, ((double[]){0.0, 0.5, 0.5, 1.0}),
        sizeof(*args->c_vec) * (ns));

    args->num_steps = 2;
    args->num_forw_sens = dims->nx + dims->nu;
    args->sens_forw = true;
    args->sens_adj = false;
    args->sens_hess = false;

    args->submodules.forward_vde->initialize_default_args(args->forward_vde_args);
    args->submodules.adjoint_vde->initialize_default_args(args->adjoint_vde_args);
    args->submodules.hess_vde->initialize_default_args(args->hess_vde_args);
}


int sim_erk_integrator_calculate_memory_size(sim_dims *dims, void *args_)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *)args_;

    int size = sizeof(sim_erk_integrator_memory);

    if (args->sens_forw) {
        size += args->submodules.forward_vde->calculate_memory_size(&sim_erk_forward_vde_dims, args->forward_vde_args);
    } 
    
    if (args->sens_adj) {
        size += args->submodules.adjoint_vde->calculate_memory_size(&sim_erk_adjoint_vde_dims, args->adjoint_vde_args);
    }

    if (args->sens_hess) {
        size += args->submodules.hess_vde->calculate_memory_size(&sim_erk_hess_vde_dims, args->hess_vde_args);
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_erk_integrator_assign_memory(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *)args_;

    sim_erk_integrator_memory *mem;

    char *c_ptr = (char *) raw_memory;

    mem = (sim_erk_integrator_memory *) c_ptr;
    c_ptr += sizeof(sim_erk_integrator_memory);

    if (args->sens_forw) {
        mem->forward_vde_mem = args->submodules.forward_vde->assign_memory(&sim_erk_forward_vde_dims, args->forward_vde_args, c_ptr);
        c_ptr += args->submodules.forward_vde->calculate_memory_size(&sim_erk_forward_vde_dims, args->forward_vde_args);
    } 

    if (args->sens_adj) {
        mem->adjoint_vde_mem = args->submodules.adjoint_vde->assign_memory(&sim_erk_adjoint_vde_dims, args->adjoint_vde_args, c_ptr);
        c_ptr += args->submodules.adjoint_vde->calculate_memory_size(&sim_erk_adjoint_vde_dims, args->adjoint_vde_args);
    }

    if (args->sens_hess) {
        mem->hess_vde_mem = args->submodules.hess_vde->assign_memory(&sim_erk_hess_vde_dims, args->hess_vde_args, c_ptr);
        c_ptr += args->submodules.hess_vde->calculate_memory_size(&sim_erk_hess_vde_dims, args->hess_vde_args);
    }

    align_char_to(8, &c_ptr);

    assert((size_t)c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    assert((char*)raw_memory + sim_erk_integrator_calculate_memory_size(dims, args_) >= c_ptr);

    return (void *)mem;
}



int sim_erk_integrator_calculate_workspace_size(sim_dims *dims, void *args_)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;
    int NF = args->num_forw_sens;

    int num_stages = args->num_stages; // number of stages
    int nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
    int nhess = (NF + 1) * NF / 2;
    uint num_steps = args->num_steps;  // number of steps

    int size = sizeof(sim_erk_integrator_workspace);

    size += (nX + nu + np) * sizeof(double); // rhs_forw_in

    if(args->sens_adj){
        size += num_steps * num_stages * nX * sizeof(double); // K_traj
        size += (num_steps + 1) * nX *sizeof(double); // out_forw_traj
    }else{
        size += num_stages * nX * sizeof(double); // K_traj
        size += nX *sizeof(double); // out_forw_traj
    }

    if (args->sens_hess && args->sens_adj){
        size += (nx + nX + nu + np) * sizeof(double); //rhs_adj_in
        size += (nx + nu + nhess) * sizeof(double); //out_adj_tmp
        size += num_stages * (nx + nu + nhess) * sizeof(double); //adj_traj
    }else if (args->sens_adj){
        size += (nx * 2 + nu + np) * sizeof(double); //rhs_adj_in
        size += (nx + nu)* sizeof(double); //out_adj_tmp
        size += num_stages * (nx + nu) * sizeof(double); //adj_traj
    }

    if (args->sens_forw) {
        size += args->submodules.forward_vde->calculate_workspace_size(&sim_erk_forward_vde_dims, args->forward_vde_args);
    } 
    
    if (args->sens_adj) {
        size += args->submodules.adjoint_vde->calculate_workspace_size(&sim_erk_adjoint_vde_dims, args->adjoint_vde_args);
    }

    if (args->sens_hess) {
        size += args->submodules.hess_vde->calculate_workspace_size(&sim_erk_hess_vde_dims, args->hess_vde_args);
    }

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}

static void *cast_workspace(sim_dims *dims, void *args_, void *raw_memory)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *) args_;

    int nx = dims->nx;
    int nu = dims->nu;
    int np = dims->np;
    int NF = args->num_forw_sens;

    int num_stages = args->num_stages; // number of stages
    int nX = nx*(1+NF); // (nx) for ODE and (NF*nx) for VDE
    int nhess = (NF + 1) * NF / 2;
    int num_steps = args->num_steps;  // number of steps

    char *c_ptr = (char *)raw_memory;

    sim_erk_integrator_workspace *workspace = (sim_erk_integrator_workspace *) c_ptr;
    c_ptr += sizeof(sim_erk_integrator_workspace);

    align_char_to(8, &c_ptr);

    assign_double(nX + nu + np, &workspace->rhs_forw_in, &c_ptr);

    if(args->sens_adj)
    {
        assign_double(num_stages*num_steps*nX, &workspace->K_traj, &c_ptr);
        assign_double((num_steps + 1)*nX, &workspace->out_forw_traj, &c_ptr);
    } else
    {
        assign_double(num_stages*nX, &workspace->K_traj, &c_ptr);
        assign_double(nX, &workspace->out_forw_traj, &c_ptr);
    }

    if (args->sens_hess && args->sens_adj)
    {
        assign_double(nx+nX+nu+np, &workspace->rhs_adj_in, &c_ptr);
        assign_double(nx+nu+nhess, &workspace->out_adj_tmp, &c_ptr);
        assign_double(num_stages*(nx+nu+nhess), &workspace->adj_traj, &c_ptr);
    } else if (args->sens_adj)
    {
        assign_double((nx*2+nu+np), &workspace->rhs_adj_in, &c_ptr);
        assign_double(nx+nu, &workspace->out_adj_tmp, &c_ptr);
        assign_double(num_stages*(nx+nu), &workspace->adj_traj, &c_ptr);
    }

    if (args->sens_forw) {
        workspace->forward_vde_work = (void *) c_ptr;
        c_ptr += args->submodules.forward_vde->calculate_workspace_size(&sim_erk_forward_vde_dims, args->forward_vde_args);
    } 

    if (args->sens_adj) {
        workspace->adjoint_vde_work = (void *) c_ptr;
        c_ptr += args->submodules.adjoint_vde->calculate_workspace_size(&sim_erk_adjoint_vde_dims, args->adjoint_vde_args);
    }

    if (args->sens_hess) {
        workspace->hess_vde_work = (void *) c_ptr;
        c_ptr += args->submodules.hess_vde->calculate_workspace_size(&sim_erk_hess_vde_dims, args->hess_vde_args);
    }

    assert((char*)raw_memory + sim_erk_integrator_calculate_workspace_size(dims, args_) >= c_ptr);

    return (void *)workspace;
}



static void compute_forward_vde(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_erk_integrator_args *args,
                         sim_erk_integrator_memory *mem,
                         sim_erk_integrator_workspace *work) 
{
    double *x = in;
    double *Sx = in + nx;
    double *Su = in + nx + nx * nx;
    double *u = in + nx + nx * (nx + nu);
    double *p = in + nx + nx * (nx + nu) + nu;

    double *x_out = out;
    double *Sx_out = out + nx;
    double *Su_out = out + nx + nx * nx;

    double *fw_in_inputs[FW_VDE_NUMIN];
    fw_in_inputs[0] = x;
    fw_in_inputs[1] = Sx;
    fw_in_inputs[2] = Su;
    fw_in_inputs[3] = u;
    fw_in_inputs[4] = p;

    bool fw_in_compute_output[FW_VDE_NUMOUT] = {true, true, true};

    double *fw_out_outputs[FW_VDE_NUMOUT];
    fw_out_outputs[0] = x_out;
    fw_out_outputs[1] = Sx_out;
    fw_out_outputs[2] = Su_out;

    external_function_in fw_vde_in;
    fw_vde_in.inputs = fw_in_inputs;
    fw_vde_in.compute_output = fw_in_compute_output;

    external_function_out fw_vde_out;
    fw_vde_out.outputs = fw_out_outputs;

    args->submodules.forward_vde->fun(&fw_vde_in, &fw_vde_out, args->forward_vde_args, mem->forward_vde_mem, work->forward_vde_work);
}



static void compute_adjoint_vde(const int nx, const int nu, const int np, double *in, double *out, 
                         sim_erk_integrator_args *args,
                         sim_erk_integrator_memory *mem,
                         sim_erk_integrator_workspace *work)
{
    double *x = in;
    double *lambdaX = in + nx;
    double *u = in + 2 * nx;
    double *p = in + 2 * nx + nu;

    double *adj_out = out;

    double *adj_in_inputs[ADJ_VDE_NUMIN];
    adj_in_inputs[0] = x;
    adj_in_inputs[1] = lambdaX;
    adj_in_inputs[2] = u;
    adj_in_inputs[3] = p;

    bool adj_in_compute_output[ADJ_VDE_NUMOUT] = {true};

    double *adj_out_outputs[ADJ_VDE_NUMOUT];
    adj_out_outputs[0] = adj_out;

    external_function_in adj_vde_in;
    adj_vde_in.inputs = adj_in_inputs;
    adj_vde_in.compute_output = adj_in_compute_output;

    external_function_out adj_vde_out;
    adj_vde_out.outputs = adj_out_outputs;

    args->submodules.adjoint_vde->fun(&adj_vde_in, &adj_vde_out, args->adjoint_vde_args, mem->adjoint_vde_mem, work->adjoint_vde_work);
}

static void compute_hess_vde(const int nx, const int nu, const int np, double *in, double *out, 
                      sim_erk_integrator_args *args,
                      sim_erk_integrator_memory *mem,
                      sim_erk_integrator_workspace *work)
{
    double *x = in;
    double *Sx = in + nx;
    double *Su = in + nx + nx * nx;
    double *lambdaX = in + nx * (1 + nx + nu);
    double *u = in + nx * (2 + nx + nu);
    double *p = in + nx * (2 + nx + nu) + nu;

    double *adj_out = out;
    double *hess_out = out + nx + nu;

    double *hess_in_inputs[HESS_VDE_NUMIN];
    hess_in_inputs[0] = x;
    hess_in_inputs[1] = Sx;
    hess_in_inputs[2] = Su;
    hess_in_inputs[3] = lambdaX;
    hess_in_inputs[4] = u;
    hess_in_inputs[5] = p;

    bool hess_in_compute_output [HESS_VDE_NUMOUT] = {true, true};

    double *hess_out_outputs[HESS_VDE_NUMOUT];
    hess_out_outputs[0] = adj_out;
    hess_out_outputs[1] = hess_out;

    external_function_in hess_vde_in;
    hess_vde_in.inputs = hess_in_inputs;
    hess_vde_in.compute_output = hess_in_compute_output;

    external_function_out hess_vde_out;
    hess_vde_out.outputs = hess_out_outputs;

    args->submodules.hess_vde->fun(&hess_vde_in, &hess_vde_out, args->hess_vde_args, mem->hess_vde_mem, work->hess_vde_work);
}



int sim_erk_integrator(sim_in *in, sim_out *out, void *args_, void *mem_, void *work_)
{
    sim_erk_integrator_args *args = (sim_erk_integrator_args *) args_;
    sim_dims dims = {
        args->num_stages,
        in->nx,
        in->nu,
        in->np
    };
    sim_erk_integrator_memory *memory = (sim_erk_integrator_memory *) mem_;
    sim_erk_integrator_workspace *workspace = (sim_erk_integrator_workspace *) cast_workspace(&dims, args, work_);

    int i, j, s, istep;
    double a = 0, b =0; // temp values of A_mat and b_vec
    int nx = in->nx;
    int nu = in->nu;
    int np = in->np;

    int NF = args->num_forw_sens;
    if (!args->sens_forw)
        NF = 0;

    int nhess = (NF + 1) * NF / 2;
    int nX = nx * (1 + NF);

    double *x = in->x;
    double *u = in->u;
    double *p = in->p;
    double *S_forw_in = in->S_forw;
    int num_steps = args->num_steps;
    double step = in->step;

    double *S_adj_in = in->S_adj;

    double *A_mat = args->A_mat;
    double *b_vec = args->b_vec;
    //    double *c_vec = args->c_vec;
    int num_stages = args->num_stages;

    double *K_traj = workspace->K_traj;
    double *forw_traj = workspace->out_forw_traj;
    double *rhs_forw_in = workspace->rhs_forw_in;

    double *adj_tmp = workspace->out_adj_tmp;
    double *adj_traj = workspace->adj_traj;
    double *rhs_adj_in = workspace->rhs_adj_in;

    double *xn = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;
    double *S_hess_out = out->S_hess;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    acados_tic(&timer);
    for (i = 0; i < nx; i++)
        forw_traj[i] = x[i];  // x0
    if (args->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            forw_traj[nx + i] = S_forw_in[i];  // sensitivities
    }

    for (i = 0; i < nu; i++)
        rhs_forw_in[nX + i] = u[i]; // controls

    for (i = 0; i < np; i++) 
        rhs_forw_in[nX + nu + i] = p[i];  // parameters

    // FORWARD SWEEP:
    for (istep = 0; istep < num_steps; istep++) {
        if (args->sens_adj) {
            K_traj = workspace->K_traj + istep * num_stages * nX;
            forw_traj = workspace->out_forw_traj + (istep + 1) * nX;
            for (i = 0; i < nX; i++)
                forw_traj[i] = forw_traj[i - nX];
        }

        for (s = 0; s < num_stages; s++) {
            for (i = 0; i < nX; i++)
                rhs_forw_in[i] = forw_traj[i];
            for (j = 0; j < s; j++){
                a = A_mat[j * num_stages + s];
                if (a!=0){
                    a *= step;
                    for (i = 0; i < nX; i++)
                        rhs_forw_in[i] += a * K_traj[j * nX + i];
                }
            }

            acados_tic(&timer_ad);
            compute_forward_vde(nx, nu, np, rhs_forw_in, K_traj+s*nX, args, memory, workspace);
            timing_ad += acados_toc(&timer_ad);
        }
        
        for (s = 0; s < num_stages; s++){
            b = step * b_vec[s];
            for (i = 0; i < nX; i++)
                forw_traj[i] += b * K_traj[s * nX + i];  // ERK step
        }
    }

    for (i = 0; i < nx; i++)
        xn[i] = forw_traj[i];
    if (args->sens_forw) {
        for (i = 0; i < nx * NF; i++)
            S_forw_out[i] = forw_traj[nx + i];
    }

    // ADJOINT SWEEP:
    if (args->sens_adj) {
        for (i = 0; i < nx; i++)
            adj_tmp[i] = S_adj_in[i];
        for (i = 0; i < nu; i++)
            adj_tmp[nx+i] = 0.0;

        int nForw = nx;
        int nAdj = nx + nu;
        if (args->sens_hess) {
            nForw = nX;
            nAdj = nx + nu + nhess;
            for (i = 0; i < nhess; i++)
                adj_tmp[nx + nu + i] = 0.0;
        }

        printf("\nnFOrw=%d nAdj=%d\n", nForw, nAdj);

        for (i = 0; i < nu; i++)
            rhs_adj_in[nForw + nx + i] = u[i];

        for (i = 0; i < np; i++)
            rhs_adj_in[nForw + nx + nu + i] = p[i];

        for (istep = num_steps - 1; istep > -1; istep--) {
            K_traj = workspace->K_traj + istep * num_stages * nX;
            forw_traj = workspace->out_forw_traj + (istep+1) * nX;
            for (s = num_stages - 1; s > -1; s--) {
                // forward variables:
                for (i = 0; i < nForw; i++)
                    rhs_adj_in[i] = forw_traj[i]; // extract x trajectory
                for (j = 0; j < s; j++)
                    a = A_mat[j * num_stages + s];
                    if (a!=0){
                        a*=step;
                        for (i = 0; i < nForw; i++)
                            rhs_adj_in[i] += a *K_traj[j * nX + i];
                    } // plus k traj
                // adjoint variables:
                b = step * b_vec[s];
                for (i = 0; i < nx; i++)
                    rhs_adj_in[nForw + i] = b * adj_tmp[i];
                for (j = s + 1; j < num_stages; j++){
                    a = A_mat[s * num_stages + j];
                    if (a!=0){
                        a *= step;
                        for (i = 0; i < nx; i++)
                            rhs_adj_in[nForw + i] += a * adj_traj[j * nAdj + i];
                    }
                }
                acados_tic(&timer_ad);
                if (args->sens_hess){
                    compute_hess_vde(nx, nu, np, rhs_adj_in, adj_traj+s*nAdj, args, memory, workspace);
                }else{
                    compute_adjoint_vde(nx, nu, np, rhs_adj_in, adj_traj+s*nAdj, args, memory, workspace);
                }
                timing_ad += acados_toc(&timer_ad);

                // printf("\nadj_traj:\n");
                // for (int ii=0;ii<num_stages*nAdj;ii++)
                //     printf("%3.1f ", adj_traj[ii]);
            }
            for (s = 0; s < num_stages; s++)
                for (i = 0; i < nAdj; i++)
                    adj_tmp[i] += adj_traj[s * nAdj + i];  // ERK step
        }
        for (i = 0; i < nx + nu; i++)
            S_adj_out[i] = adj_tmp[i];
        if (args->sens_hess) {
            for (i = 0; i < nhess; i++)
                S_hess_out[i] = adj_tmp[nx + nu + i];
        }
    }
    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;
    return 0;  // success
}
