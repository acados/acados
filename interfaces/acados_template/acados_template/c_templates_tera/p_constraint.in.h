#ifndef {{ con_p_name }}_P_CONSTRAINT
#define {{ con_p_name }}_P_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if dims.npd > 0 %}
// implicit ODE
int {{ con_p_name }}_p_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ con_p_name }}_p_constraint_work(int *, int *, int *, int *);
const int *{{ con_p_name }}_p_constraint_sparsity_in(int);
const int *{{ con_p_name }}_p_constraint_sparsity_out(int);
int {{ con_p_name }}_p_constraint_n_in();
int {{ con_p_name }}_p_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ con_p_name }}_P_CONSTRAINT
