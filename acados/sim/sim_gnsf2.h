typedef struct
{
    int num_stages;
    int nx;
    int nu;
    int nz;
    int nx1;
    int nx2;
    int num_steps;
    int n_out;
    int n_in;
    
} gnsf2_dims;


typedef struct {
    double *ff;
    double *x0_1;
    double *u_0;
} gnsf_res_in;

typedef struct {
    double *x;
    double *u;
    double *S_forw;  // forward seed
    double *S_adj;   // backward seed

} gnsf_in;

typedef struct {
    struct blasfeo_dmat KKf;
    struct blasfeo_dmat KKx;
    struct blasfeo_dmat KKu;

    struct blasfeo_dmat ZZf;
    struct blasfeo_dmat ZZx;
    struct blasfeo_dmat ZZu;

    struct blasfeo_dmat ALO;
    struct blasfeo_dmat M2inv;
    struct blasfeo_dmat dK2_dx2;

    double* A_dt;
    double* b_dt;
    double* c;
    double dt;

    // external functions
    external_function_generic *res_inc_Jff;
    external_function_generic *jac_res_ffx1u;
    external_function_generic *f_LO_inc_J_x1k1uz;
    
} gnsf_fixed;

typedef struct
{
	/* external functions */
	// nonlinearity functions
	external_function_generic *Phi_inc_dy;
	// Linear output function
	external_function_generic *f_LO_inc_J_x1k1uz;

    // model defining matrices
    double *A;
    double *B;
    double *C;
    double *E;
    double *L_x;
    double *L_xdot;
    double *L_z;
    double *L_u;
    double *ALO;
} gnsf2_model;