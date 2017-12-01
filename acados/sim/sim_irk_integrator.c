// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/mem.h"

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"

int sim_irk_opts_calculate_size(sim_dims *dims)
{

    int size = sizeof(sim_rk_opts);

    int ns = dims->num_stages;
    size += ns * ns * sizeof(double);  // A_mat
    size += ns * sizeof(double);  // b_vec
    size += ns * sizeof(double);  // c_vec

    make_int_multiple_of(8, &size);
    size += 1 * 8;

    return size;
}



void *sim_irk_assign_opts(sim_dims *dims, void *raw_workspaceory)
{
    char *c_ptr = (char *) raw_workspaceory;

    sim_rk_opts *opts = (sim_rk_opts *) c_ptr;
    c_ptr += sizeof(sim_rk_opts);

    int ns = dims->num_stages;
    opts->num_stages = ns;

    align_char_to(8, &c_ptr);

    assign_double(ns*ns, &opts->A_mat, &c_ptr);
    assign_double(ns, &opts->b_vec, &c_ptr);
    assign_double(ns, &opts->c_vec, &c_ptr);

    assert((char*)raw_workspaceory + sim_irk_opts_calculate_size(dims) >= c_ptr);

    opts->newton_iter = 3;
    opts->scheme = NULL;

    return (void *)opts;
}


void sim_irk_initialize_default_args(sim_dims *dims, void *opts_)
{
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    int ns = opts->num_stages;

    assert(opts->num_stages == 3 && "only number of stages = 3 implemented!");

    memcpy(opts->A_mat,
        ((real_t[]){0.1389, 0.3003, 0.2680 ,
        -0.0360, 0.2222, 0.4804,
        0.0098, -0.0225, 0.1389}),
        sizeof(*opts->A_mat) * (ns * ns));
    memcpy(opts->b_vec, ((real_t[]){0.2778, 0.4444, 0.2778}),
        sizeof(*opts->b_vec) * (ns));
    memcpy(opts->c_vec, ((real_t[]){0.1127, 0.5, 0.8873}),
        sizeof(*opts->c_vec) * (ns));

    opts->num_steps = 2;
    opts->num_forw_sens = dims->nx + dims->nu;
    opts->sens_forw = true;
    opts->sens_adj = false;
    opts->sens_hess = false;
}


int sim_irk_calculate_memory_size(sim_dims *dims, void *opts_)
{
    return 0;
}

void *sim_irk_assign_memory(sim_dims *dims, void *opts_, void *raw_workspaceory)
{
    return NULL;
}

int sim_irk_calculate_workspace_size(sim_dims *dims, void *opts_)
{
    sim_rk_opts *opts = (sim_rk_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    
    int ns = opts->num_stages; // number of stages
    int steps = opts->num_steps;

    int size = sizeof(sim_irk_workspace);

    size += 5*sizeof(struct d_strmat); // JG, JGf, JKf, JFK, S_forw
    size += 4*sizeof(struct d_strvec); // rG, K, xt, xn
    size += 2*sizeof(struct d_strvec); // lambda,lambdaK

    size += 2*sizeof(struct d_strvec *); // **xn_traj, **K_traj;
    size += 2*steps*sizeof(struct d_strvec);  // number of strvec pointers

    size += d_size_strmat(nx*ns, nx*ns); // JG
    size += 2*d_size_strmat(nx*ns, nx+nu); // JGf, JKf
    size += d_size_strmat(nx, nx*ns); // JFK
    size += d_size_strmat(nx, nx+nu); // S_forw

    size += 2*d_size_strvec(nx*ns); // rG, K
    size += 2*d_size_strvec(nx); // xt, x

    size += d_size_strvec(nx+nu); // lambda
    size += d_size_strvec(nx*ns); // lambdaK

    size += steps * d_size_strvec(nx); // for xn_traj
    size += steps * d_size_strvec(nx*ns); // for K_traj

    size += nx * sizeof(double); //  rGt
    size += nx * (2*nx+nu) * sizeof(double); // jac_out
    size += nx * nx * sizeof(double); // Jt
    size += (2*nx + nu) * sizeof(double); // ode_args
    size += (nx+nu) * sizeof(double); // S_adj_w

    size += nx *ns * sizeof(int); // ipiv

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

void *sim_irk_cast_workspace(sim_dims *dims, void *opts_, void *raw_memory)
{

    sim_rk_opts *opts = (sim_rk_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
        
    int ns = opts->num_stages; // number of stages
    int steps = opts->num_steps;

    char *c_ptr = (char *)raw_memory;

    sim_irk_workspace *workspace = (sim_irk_workspace *) c_ptr;
    c_ptr += sizeof(sim_irk_workspace);

    workspace->xn_traj = (struct d_strvec **)c_ptr;
    c_ptr += steps* sizeof(struct d_strvec *); 
    
    workspace->K_traj = (struct d_strvec **)c_ptr;
    c_ptr += steps* sizeof(struct d_strvec *); 

    workspace->JG = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat); 

    workspace->JGf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    workspace->JKf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    workspace->JFK = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    workspace->S_forw = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    workspace->rG = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    workspace->K = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    workspace->xt = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    workspace->xn = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    workspace->lambda = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    workspace->lambdaK = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    for (int i=0;i<steps;i++){
        workspace->xn_traj[i] = (struct d_strvec *)c_ptr;
        c_ptr += sizeof(struct d_strvec);

        workspace->K_traj[i] = (struct d_strvec *)c_ptr;
        c_ptr += sizeof(struct d_strvec);
    }
    
    align_char_to(64, &c_ptr);

    assign_strmat(nx*ns, nx*ns, workspace->JG, &c_ptr);

    assign_strmat(nx*ns, nx+nu, workspace->JGf, &c_ptr);

    assign_strmat(nx*ns, nx+nu, workspace->JKf, &c_ptr);

    assign_strmat(nx, nx*ns, workspace->JFK, &c_ptr);

    assign_strmat(nx, nx+nu, workspace->S_forw, &c_ptr);

    assign_strvec(nx*ns, workspace->rG, &c_ptr);

    assign_strvec(nx*ns, workspace->K, &c_ptr);

    assign_strvec(nx, workspace->xt, &c_ptr);

    assign_strvec(nx, workspace->xn, &c_ptr);

    assign_strvec(nx+nu, workspace->lambda, &c_ptr);

    assign_strvec(nx*ns, workspace->lambdaK, &c_ptr);

    for (int i=0;i<steps;i++){
        assign_strvec(nx, workspace->xn_traj[i], &c_ptr);
        assign_strvec(nx*ns, workspace->K_traj[i], &c_ptr);
    }

    assign_double(nx, &workspace->rGt, &c_ptr);

    assign_double(nx * (2*nx+nu), &workspace->jac_out, &c_ptr);

    assign_double(nx * nx, &workspace->Jt, &c_ptr);

    assign_double(2*nx + nu, &workspace->ode_args, &c_ptr);

    assign_double(nx + nu, &workspace->S_adj_w, &c_ptr);

    assign_int(nx * ns , &workspace->ipiv, &c_ptr);

    assert((char*)raw_memory + sim_irk_calculate_workspace_size(dims, opts_) >= c_ptr);

    return (void *)workspace;
}

void *sim_irk_create_memory(sim_dims *dims, void *opts_)
{

    int bytes = sim_irk_calculate_memory_size(dims, opts_);
    void *ptr = malloc(bytes);
    sim_irk_memory *memory = sim_irk_assign_memory(dims, opts_, ptr);

    return (void *)memory;
}


int sim_irk(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_){
    
    sim_rk_opts *opts = (sim_rk_opts *) opts_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;

    sim_dims dims = {
        opts->num_stages,
        in->nx,
        in->nu
    };
    sim_irk_workspace *workspace = (sim_irk_workspace *) sim_irk_cast_workspace(&dims, opts, work_);

    int ii, jj, iter, kk, ss;
    double a,b;

    int nx = in->nx;
	int nu = in->nu;
    int num_steps = opts->num_steps;
    double *x = in->x;
	double *u = in->u;
    double step = in->step;
    double *S_forw_in = in->S_forw;
    
    int num_stages = opts->num_stages;
    int newton_iter = opts->newton_iter;
    double *A_mat = opts->A_mat;
    double *b_vec = opts->b_vec;
 
	double *rGt = workspace->rGt; 
	double *jac_out = workspace->jac_out;
    double *Jt = workspace->Jt;
    double *ode_args = workspace->ode_args;
    int *ipiv = workspace->ipiv;
	struct d_strmat *JG = workspace->JG;
	struct d_strvec *rG = workspace->rG;
    struct d_strvec *K = workspace->K;
    struct d_strmat *JGf = workspace->JGf;
    struct d_strmat *JKf = workspace->JKf;
    struct d_strvec *xt = workspace->xt;
    struct d_strvec *xn = workspace->xn;
    struct d_strmat *JFK = workspace->JFK;
    struct d_strmat *S_forw = workspace->S_forw;

    // for adjoint
    struct d_strvec *lambda = workspace->lambda;
    struct d_strvec *lambdaK = workspace->lambdaK;
    struct d_strvec **xn_traj = workspace->xn_traj;
    struct d_strvec **K_traj = workspace->K_traj;
    double *S_adj_in = workspace->S_adj_w;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;

    acados_timer timer, timer_ad;
    double timing_ad = 0.0;

    // initialize
    dgese_libstr(nx*num_stages, nx*num_stages, 0.0, JG, 0, 0);
    dgese_libstr(nx*num_stages, nx+nu, 0.0, JGf, 0, 0);
    dgese_libstr(nx*num_stages, nx+nu, 0.0, JKf, 0, 0);
    for (ii=0;ii<num_stages;ii++){
        b = step * b_vec[ii];
        ddiare_libstr(nx, b, JFK, 0, ii*nx);
    }
    d_cvt_mat2strmat(nx, nx+nu, S_forw_in, nx, S_forw, 0, 0);

    dvecse_libstr(nx*num_stages, 0.0, rG, 0);
    dvecse_libstr(nx*num_stages, 0.0, K, 0);
    d_cvt_vec2strvec(nx, x, xn, 0);

    dvecse_libstr(nx*num_stages, 0.0, lambdaK, 0);

    for(kk=0;kk<nx;kk++)
        S_adj_in[kk] = in->S_adj[kk];
    for(kk=0;kk<nu;kk++)
        S_adj_in[nx+kk] = 0.0;
    d_cvt_vec2strvec(nx+nu, S_adj_in, lambda, 0);
    
    for (kk=0;kk<2*nx;kk++) //initialize x,xdot with zeros
        ode_args[kk] = 0.0;
    for (kk=0;kk<nu;kk++) //set controls
        ode_args[2*nx+kk] = u[kk];

    for (kk=0;kk<nx;kk++){
        rGt[kk] = 0.0;
    }
    for (kk=0;kk<nx*(2*nx+nu);kk++)
        jac_out[kk] = 0.0;
    for (kk=0;kk<nx*nx;kk++)
        Jt[kk] = 0.0;
    for (kk=0;kk<num_stages*nx;kk++)
        ipiv[kk] = kk; //ipiv?
        
    // start the loop
    acados_tic(&timer);
    for(ss=0; ss<num_steps; ss++){

        //  obtain Kn
        for(iter=0; iter<newton_iter; iter++){

            if (opts->sens_adj){
                dveccp_libstr(nx, xn, 0, xn_traj[ss], 0);
            }
           
            for(ii=0; ii<num_stages; ii++){ // ii-th row of tableau   
                
                // take x(n); copy a strvec into a strvec
                dveccp_libstr(nx, xn, 0, xt, 0);

                for(jj=0; jj<num_stages; jj++){ // jj-th col of tableau
                    a = A_mat[ii+num_stages*jj];
					if(a!=0){           // xt = xt + T_int * a[i,j]*K_j
                        a *= step; 
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0); 
					}
                }
                // put xn+sum kj into first nx elements of ode_arg               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);
                
                // put ki into the next nx elements of ode_args
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);      

                // compute the residual of implicit ode at time t_ii, store value in rGt 
                in->eval_impl_res(nx, nu, ode_args, rGt, in->impl_ode);
                // fill in elements of rG  - store values rGt on (ii*nx)th position of rG
                d_cvt_vec2strvec(nx, rGt, rG, ii*nx);
                
                // compute the jacobian of implicit ode
                acados_tic(&timer_ad);
                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                timing_ad += acados_toc(&timer_ad);
                // compute the blocks of JG
                for (jj=0; jj<num_stages; jj++){ //compute the block (ii,jj)th block = Jt
                    a = A_mat[ii + num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0; kk<nx*nx; kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    // fill in the ii-th, jj-th block of JG
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);
                } // end jj

            } // end ii
            
            //DGETRF computes an LU factorization of a general M-by-N matrix A
            //using partial pivoting with row interchanges. // and store it in JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);
            // permute also the r.h.s
            dvecpe_libstr(nx*num_stages, ipiv, rG, 0);

            // solve JG * y = rG, JG on the (l)eft, (l)ower-trian, (n)o-trans
            //                    (u)nit trian
            dtrsv_lnu_libstr(nx*num_stages, JG, 0, 0, rG, 0, rG, 0);
            
            // solve JG * x = rG, JG on the (l)eft, (u)pper-trian, (n)o-trans
            //                    (n)o unit trian , and store x in rG
            dtrsv_unn_libstr(nx*num_stages, JG, 0, 0, rG, 0, rG, 0);
            // scale and add a generic strmat into a generic strmat // K = K - rG, where rG is DeltaK
            daxpy_libstr(nx*num_stages, -1.0, rG, 0, K, 0, K, 0);
        }// end iter

        if (opts->sens_adj){
            dveccp_libstr(nx*num_stages, K, 0, K_traj[ss], 0);
        }
        
        // evaluate forward sensitivities
        if (opts->sens_forw){

            // evaluate JG(xn,Kn)         
            for(ii=0; ii<num_stages; ii++){  

                dveccp_libstr(nx, xn, 0, xt, 0);
                //compute xt = final x;
                for(jj=0; jj<num_stages; jj++){ 
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){
                        a *= step;
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0);
                    }
                }               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);                                           
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);  
                
                acados_tic(&timer_ad);
                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_u(nx, nu, ode_args, jac_out+2*nx*nx, in->impl_jac_u);
                d_cvt_mat2strmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                timing_ad += acados_toc(&timer_ad);
                // maybe this part can also be written using blasfeo
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);            
                } // end jj     
            } // end ii

            // factorize JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);

            // obtain JKf
            dgemm_nn_libstr(nx*num_stages, nx, nx, 1.0, JGf, 0, 0, S_forw, 0, 0, 0.0, JKf, 0, 0, JKf, 0, 0);
            dgemm_nn_libstr(nx*num_stages, nu, nx, 1.0, JGf, 0, 0, S_forw, 0, nx, 1.0, JGf, 0, nx, JKf, 0, nx);  
            drowpe_libstr(nx*num_stages, ipiv, JKf);
            dtrsm_llnu_libstr(nx*num_stages, nx+nu, 1.0, JG, 0, 0, JKf, 0, 0, JKf, 0, 0); 
            dtrsm_lunn_libstr(nx*num_stages, nx+nu, 1.0, JG, 0, 0, JKf, 0, 0, JKf, 0, 0);

            // update forward sensitivity
            dgemm_nn_libstr(nx, nx+nu, nx*num_stages, -1.0, JFK, 0, 0, JKf, 0, 0, 1.0, S_forw, 0, 0, S_forw, 0, 0);
        }

        // obtain x(n+1)
        for(ii=0;ii<num_stages;ii++){
            b = step * b_vec[ii];
            daxpy_libstr(nx, b, K, ii*nx, xn, 0, xn, 0);          
        }
            
    }// end int step ss

    // evaluate backwards
    if (opts->sens_adj){ 

        for(ss=num_steps-1;ss>-1;ss--){

            dveccp_libstr(nx, xn_traj[ss], 0, xt, 0);

            dveccp_libstr(num_stages*nx, K_traj[ss], 0, K, 0);

            for(ii=0; ii<num_stages; ii++){  

                for(jj=0; jj<num_stages; jj++){ 
                    a = A_mat[ii+num_stages*jj];
                    if(a!=0){
                        a *= step;
                        daxpy_libstr(nx, a, K, jj*nx, xt, 0, xt, 0);
                    }
                }               
                d_cvt_strvec2vec(nx, xt, 0, ode_args);                                           
                d_cvt_strvec2vec(nx, K, ii*nx, ode_args+nx);  

                acados_tic(&timer_ad);
                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_u(nx, nu, ode_args, jac_out+2*nx*nx, in->impl_jac_u);
                d_cvt_mat2strmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                timing_ad += acados_toc(&timer_ad);
                for (jj=0;jj<num_stages;jj++){
                    a = A_mat[ii+num_stages*jj];
                    if (a!=0){
                        a *= step;
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] = a* jac_out[kk];                          
                    }
                    if(jj==ii){
                        for (kk=0;kk<nx*nx;kk++)
                            Jt[kk] += jac_out[nx*nx+kk];
                    }
                    d_cvt_mat2strmat(nx, nx, Jt, nx, JG, ii*nx, jj*nx);            
                } // end jj     
            } // end ii

            // factorize JG
            dgetrf_libstr(nx*num_stages, nx*num_stages, JG, 0, 0, JG, 0, 0, ipiv);

            dgemv_t_libstr(nx, nx*num_stages, -1.0, JFK, 0, 0, lambda, 0, 0.0, lambdaK, 0, lambdaK, 0);

            dtrsv_utn_libstr(nx*num_stages, JG, 0, 0, lambdaK, 0, lambdaK, 0);

            dtrsv_ltu_libstr(nx*num_stages, JG, 0, 0, lambdaK, 0, lambdaK, 0);

            dvecpei_libstr(nx*num_stages, ipiv, lambdaK, 0);

            // ???
            // dgemv_t_libstr(nx*num_stages, nx, 1.0, JGf, 0, 0, lambdaK, 0, 1.0, lambda, 0, lambda, 0);
            // dgemv_t_libstr(nx*num_stages, nu, 1.0, JGf, 0, nx, lambdaK, 0, 1.0, lambda, nx, lambda, nx);
            dgemv_t_libstr(nx*num_stages, nx+nu, 1.0, JGf, 0, 0, lambdaK, 0, 1.0, lambda, 0, lambda, 0);
        }
    }

    // extract output

    d_cvt_strvec2vec(nx, xn, 0, x_out);

    d_cvt_strmat2mat(nx, nx+nu, S_forw, 0, 0, S_forw_out, nx);

    d_cvt_strvec2vec(nx+nu, lambda, 0, S_adj_out);

    out->info->CPUtime = acados_toc(&timer);
    out->info->LAtime = 0.0;
    out->info->ADtime = timing_ad;

    return 0;
}