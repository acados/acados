#ifndef {{ ra.con_h_name.upper() }}_H_CONSTRAINT
#define {{ ra.con_h_name.upper() }}_H_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if ra.dims.nh > 0: %}
// implicit ODE
int {{ ra.con_h_name }}_h_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.con_h_name }}_h_constraint_work(int *, int *, int *, int *);
const int *{{ ra.con_h_name }}_h_constraint_sparsity_in(int);
const int *{{ ra.con_h_name }}_h_constraint_sparsity_out(int);
int {{ ra.con_h_name }}_h_constraint_n_in();
int {{ ra.con_h_name }}_h_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ra.con_h_name.upper() }}_H_CONSTRAINT
