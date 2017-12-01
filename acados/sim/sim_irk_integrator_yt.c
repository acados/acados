// standard
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
// acados
#include "acados/utils/print.h"
#include "acados/utils/mem.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_irk_integrator_yt.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"


int irk_calculate_memory_size(sim_dims *dims, void *opts_)
{
    sim_RK_opts *opts = (sim_RK_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
    
    int ns = opts->num_stages; // number of stages
    int steps = opts->num_steps;

    int size = sizeof(sim_irk_memory);

    size += 5*sizeof(struct d_strmat); // JG, JGf, JKf, JFK
    size += 4*sizeof(struct d_strvec); // rG, K, xt, xn
    size += 2*sizeof(struct d_strvec); // lambda,lambdaK

    size += 2*sizeof(struct d_strvec *); // **xn_traj, **K_traj;
    size += 2*steps*sizeof(struct d_strvec);  // number of strvec pointers

    size += d_size_strmat(nx*ns, nx*ns); // JG
    size += 2*d_size_strmat(nx*ns, nx+nu); // JGf, JKf
    size += d_size_strmat(nx, nx*ns); // JFK
    size += d_size_strmat(nx, nx+nu); // S_forw

    size += 4*d_size_strvec(nx*ns); // rG, K
    size += 2*d_size_strvec(nx); // xt, x

    size += d_size_strvec(nx+nu); // lambda
    size += d_size_strvec(nx*ns); // lambdaK

    size += steps * d_size_strvec(nx); // for xn_traj
    size += steps * d_size_strvec(nx*ns); // for K_traj

    size += nx * sizeof(double); //  rGt
    size += nx * (2*nx+nu) * sizeof(double); // jac_out
    size += nx * nx * sizeof(double); // Jt
    size += (2*nx + nu) * sizeof(double); // ode_args
    size += (nx+nu) * sizeof(double);

    size += nx *ns * sizeof(int); // ipiv

    make_int_multiple_of(64, &size);
    size += 1 * 64;

    return size;
}

char *assign_irk_memory(sim_dims *dims, void *opts_, void *raw_memory)
{

    sim_RK_opts *opts = (sim_RK_opts *) opts_;

    int nx = dims->nx;
    int nu = dims->nu;
        
    int ns = opts->num_stages; // number of stages
    int steps = opts->num_steps;

    char *c_ptr = (char *)raw_memory;

    sim_irk_memory *mem = (sim_irk_memory *) c_ptr;
    c_ptr += sizeof(sim_irk_memory);


    (*memory)->xn_traj = (struct d_strvec **)c_ptr;
    c_ptr += steps* sizeof(struct d_strvec *); 
    
    (*memory)->K_traj = (struct d_strvec **)c_ptr;
    c_ptr += steps* sizeof(struct d_strvec *); 

    (*memory)->JG = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat); 

    (*memory)->JGf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->JKf = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->JFK = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->S_forw = (struct d_strmat *)c_ptr;
    c_ptr += sizeof(struct d_strmat);

    (*memory)->rG = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->K = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->xt = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

     (*memory)->xn = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

     (*memory)->lambda = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);

    (*memory)->lambdaK = (struct d_strvec *)c_ptr;
    c_ptr += sizeof(struct d_strvec);


    //
    for (int i=0;i<steps;i++){
        (*memory)->xn_traj[i] = (struct d_strvec *)c_ptr;
        c_ptr += sizeof(struct d_strvec);

        (*memory)->K_traj[i] = (struct d_strvec *)c_ptr;
        c_ptr += sizeof(struct d_strvec);
    }
    
    // struct d_strmat *JG = (*memory)->JG; 
    // struct d_strmat *JGf = (*memory)->JGf;
    // struct d_strmat *JKf = (*memory)->JKf;
    // struct d_strmat *JFK = (*memory)->JFK;
    // struct d_strmat *S_forw = (*memory)->S_forw;

    // struct d_strvec *rG = (*memory)->rG;
    // struct d_strvec *K = (*memory)->K;
    // struct d_strvec *xt = (*memory)->xt;
    // struct d_strvec *xn = (*memory)->xn;

    // struct d_strvec *lambda = (*memory)->lambda;
    // struct d_strvec *lambdaK = (*memory)->lambdaK;

    // struct d_strvec **xn_traj = (*memory)->xn_traj;
    // struct d_strvec **K_traj = (*memory)->K_traj;

    // align memory to typical cache line size
    // size_t s_ptr = (size_t)c_ptr;
    // s_ptr = (s_ptr + 63) / 64 * 64;
    // c_ptr = (char *)s_ptr;
    align_char_to(64, &c_ptr);

    // d_create_strmat(nx*ns, nx*ns, JG, c_ptr);
    // c_ptr += d_size_strmat(nx*ns, nx*ns);

    assign_strmat(nx*ns, nx*ns, &mem->JG, &c_ptr);

    // d_create_strmat(nx*ns, nx+nu, JGf, c_ptr);
    // c_ptr += d_size_strmat(nx*ns, nx+nu);

    assign_strmat(nx*ns, nx+nu, &mem->JGf, &c_ptr);

    // d_create_strmat(nx*ns, nx+nu, JKf, c_ptr);
    // c_ptr += d_size_strmat(nx*ns, nx+nu);

    assign_strmat(nx*ns, nx+nu, &mem->JKf, &c_ptr);

    // d_create_strmat(nx, nx*ns, JFK, c_ptr);
    // c_ptr += d_size_strmat(nx, nx*ns);

    assign_strmat(nx, nx*ns, &mem->JFK, &c_ptr);

    // d_create_strmat(nx, nx+nu, S_forw, c_ptr);
    // c_ptr += d_size_strmat(nx, nx+nu);

    assign_strmat(nx, nx+nu, &mem->S_forw, &c_ptr);

    // d_create_strvec(nx*ns, rG, c_ptr);
    // c_ptr += d_size_strvec(nx*ns);

    assign_strvec(nx*ns, &mem->rG, &c_ptr);

    // d_create_strvec(nx*ns, K, c_ptr);
    // c_ptr += d_size_strvec(nx*ns);

    assign_strvec(nx*ns, &mem->K, &c_ptr);

    // d_create_strvec(nx, xt, c_ptr);
    // c_ptr += d_size_strvec(nx);

    assign_strvec(nx, &mem->xt, &c_ptr);

    // d_create_strvec(nx, xn, c_ptr);
    // c_ptr += d_size_strvec(nx);

    assign_strvec(nx, &mem->xn, &c_ptr);

    // d_create_strvec(nx+nu, lambda, c_ptr);
    // c_ptr += d_size_strvec(nx);

    assign_strvec(nx+nu, &mem->lambda, &c_ptr);

    // d_create_strvec(nx*ns, lambdaK, c_ptr);
    // c_ptr += d_size_strvec(nx*ns);

    assign_strvec(nx*ns, &mem->lambdaK, &c_ptr);

    //
    for (int i=0;i<steps;i++){
        // d_create_strvec(nx, xn_traj[i], c_ptr);
        // c_ptr += d_size_strvec(nx);

        // d_create_strvec(nx*ns, K_traj[i], c_ptr);
        // c_ptr += d_size_strvec(nx*ns);

        assign_strvec(nx, &mem->xn_traj[i], &c_ptr);
        assign_strvec(nx*ns, &mem->K_traj[i], &c_ptr);
    }

    //
    // (*memory)->rGt = (double *)c_ptr;
    // c_ptr += nx * sizeof(double);
    assign_double(nx, &mem->rGt, &c_ptr);

    // (*memory)->jac_out = (double *)c_ptr;
    // c_ptr += nx * (2*nx+nu) * sizeof(double);

    assign_double(nx * (2*nx+nu), &mem->jac_out, &c_ptr);

    // (*memory)->Jt = (double *)c_ptr;
    // c_ptr += nx * nx * sizeof(double);

    assign_double(nx * nx, &mem->Jt, &c_ptr);

    // (*memory)->ode_args = (double *)c_ptr;
    // c_ptr += (2*nx + nu) * sizeof(double);

    assign_double(2*nx + nu, &mem->ode_args, &c_ptr);

    // (*memory)->S_adj_w = (double *)c_ptr;
    // c_ptr += (nx + nu) * sizeof(double);

    assign_double(nx + nu, &mem->S_adj_w, &c_ptr);

    // (*memory)->ipiv = (int *)c_ptr;
    // c_ptr += nx * ns * sizeof(int);

    assign_int(nx * ns , &mem->ipiv, &c_ptr);

    assert((char*)raw_memory + irk_calculate_memory_size(dims, opts_) >= c_ptr);

    return (void *)mem;
}

void *sim_irk_create_memory(sim_dims *dims, void *opts_)
{

    int bytes = irk_calculate_memory_size(dims, opts_);
    void *ptr = malloc(bytes);
    sim_irk_memory *memory = assign_irk_memory(dims, opts_, ptr);

    return (void *)memory;
}


int sim_irk_yt(sim_in *in, sim_out *out, void *opts_, void *mem_, void *work_){
    
    sim_RK_opts *opts = (sim_RK_opts *) opts_;
    sim_irk_memory *mem = (sim_irk_memory *) mem_;
    sim_irk_workspace *work = (sim_irk_workspace *) work_;

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
 
	double *rGt = mem->rGt;  // residual of one stage of the impl. ODE
	double *jac_out = mem->jac_out;
    double *Jt = mem->Jt;
    double *ode_args = mem->ode_args;
    int *ipiv = mem->ipiv;
	struct d_strmat *JG = mem->JG;
	struct d_strvec *rG = mem->rG; // residual of the of the IRK equations G
    struct d_strvec *K = mem->K;
    struct d_strmat *JGf = mem->JGf;
    struct d_strmat *JKf = mem->JKf;
    struct d_strvec *xt = mem->xt;
    struct d_strvec *xn = mem->xn;
    struct d_strmat *JFK = mem->JFK;
    struct d_strmat *S_forw = mem->S_forw;

    // for adjoint
    struct d_strvec *lambda = mem->lambda;
    struct d_strvec *lambdaK = mem->lambdaK;
    struct d_strvec **xn_traj = mem->xn_traj;
    struct d_strvec **K_traj = mem->K_traj;
    double *S_adj_in = mem->S_adj_w;

    double *x_out = out->xn;
    double *S_forw_out = out->S_forw;
    double *S_adj_out = out->S_adj;

    acados_timer timer;

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

            if (in->sens_adj){
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
                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
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

        if (in->sens_adj){
            dveccp_libstr(nx*num_stages, K, 0, K_traj[ss], 0);
        }
        
        // evaluate forward sensitivities
        if (in->sens_forw){

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

                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_u(nx, nu, ode_args, jac_out+2*nx*nx, in->impl_jac_u);
                d_cvt_mat2strmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
                
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
    if (in->sens_adj){ 

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

                in->eval_impl_jac_x(nx, nu, ode_args, jac_out, in->impl_jac_x);
                d_cvt_mat2strmat(nx, nx, jac_out, nx, JGf, ii*nx, 0);

                in->eval_impl_jac_u(nx, nu, ode_args, jac_out+2*nx*nx, in->impl_jac_u);
                d_cvt_mat2strmat(nx, nu, jac_out+2*nx*nx, nx, JGf, ii*nx, nx);

                in->eval_impl_jac_xdot(nx, nu, ode_args, jac_out+nx*nx, in->impl_jac_xdot);
            
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

    out->info->CPUtime = acados_toc(&timer);

    d_cvt_strvec2vec(nx, xn, 0, x_out);

    d_cvt_strmat2mat(nx, nx+nu, S_forw, 0, 0, S_forw_out, nx);

    d_cvt_strvec2vec(nx+nu, lambda, 0, S_adj_out);

    return 0;
}