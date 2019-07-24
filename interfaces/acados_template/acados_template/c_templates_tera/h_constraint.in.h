#ifndef {{ con_h.name }}_H_CONSTRAINT
#define {{ con_h.name }}_H_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if dims.nh > 0 %}
// implicit ODE
int {{ con_h.name }}_h_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ con_h.name }}_h_constraint_work(int *, int *, int *, int *);
const int *{{ con_h.name }}_h_constraint_sparsity_in(int);
const int *{{ con_h.name }}_h_constraint_sparsity_out(int);
int {{ con_h.name }}_h_constraint_n_in();
int {{ con_h.name }}_h_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ con_h.name }}_H_CONSTRAINT
