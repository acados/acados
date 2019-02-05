#ifndef {{ ra.con_p_name.upper() }}_P_CONSTRAINT
#define {{ ra.con_p_name.upper() }}_P_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if ra.dims.npd > 0: %}
// implicit ODE
int {{ ra.con_p_name }}_p_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ra.con_p_name }}_p_constraint_work(int *, int *, int *, int *);
const int *{{ ra.con_p_name }}_p_constraint_sparsity_in(int);
const int *{{ ra.con_p_name }}_p_constraint_sparsity_out(int);
int {{ ra.con_p_name }}_p_constraint_n_in();
int {{ ra.con_p_name }}_p_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ra.con_p_name.upper() }}_P_CONSTRAINT
