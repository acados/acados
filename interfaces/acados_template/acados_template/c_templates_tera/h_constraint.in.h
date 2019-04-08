#ifndef {{ ocp.con_h_name }}_H_CONSTRAINT
#define {{ ocp.con_h_name }}_H_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if ocp.dims.nh > 0 %}
// implicit ODE
int {{ ocp.con_h_name }}_h_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.con_h_name }}_h_constraint_work(int *, int *, int *, int *);
const int *{{ ocp.con_h_name }}_h_constraint_sparsity_in(int);
const int *{{ ocp.con_h_name }}_h_constraint_sparsity_out(int);
int {{ ocp.con_h_name }}_h_constraint_n_in();
int {{ ocp.con_h_name }}_h_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ocp.con_h_name }}_H_CONSTRAINT
