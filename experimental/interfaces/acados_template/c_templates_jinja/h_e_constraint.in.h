#ifndef {{ ocp.model.name }}_H_E_CONSTRAINT
#define {{ ocp.model.name }}_H_E_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if ocp.dims.nh_e > 0 %}
int {{ ocp.model.name }}_h_e_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ ocp.model.name }}_h_e_constraint_work(int *, int *, int *, int *);
const int *{{ ocp.model.name }}_h_e_constraint_sparsity_in(int);
const int *{{ ocp.model.name }}_h_e_constraint_sparsity_out(int);
int {{ ocp.model.name }}_h_e_constraint_n_in();
int {{ ocp.model.name }}_h_e_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ ocp.model.name }}_H_E_CONSTRAINT
