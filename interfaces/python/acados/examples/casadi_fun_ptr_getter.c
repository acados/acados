int expl_vde_for(const double** arg, double** res, int* iw, double* w, void* mem);
int expl_vde_for_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w);
const int* expl_vde_for_sparsity_in(int i);
const int* expl_vde_for_sparsity_out(int i);
int expl_vde_for_n_in(void);
int expl_vde_for_n_out(void);

void *get_fun_fun() { return &expl_vde_for; }
void *get_fun_work() { return &expl_vde_for_work; }
void *get_fun_sparsity_in() { return &expl_vde_for_sparsity_in; }
void *get_fun_sparsity_out() { return &expl_vde_for_sparsity_out; }
void *get_fun_n_in() { return &expl_vde_for_n_in; }
void *get_fun_n_out() { return &expl_vde_for_n_out; }
