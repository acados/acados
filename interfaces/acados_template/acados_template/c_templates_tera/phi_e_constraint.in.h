#ifndef {{ con_h_e.name }}_PHI_E_CONSTRAINT
#define {{ con_h_e.name }}_PHI_E_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if dims.nh_e > 0 %}
int {{ con_h_e.name }}_h_e_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ con_h_e.name }}_h_e_constraint_work(int *, int *, int *, int *);
const int *{{ con_h_e.name }}_h_e_constraint_sparsity_in(int);
const int *{{ con_h_e.name }}_h_e_constraint_sparsity_out(int);
int {{ con_h_e.name }}_h_e_constraint_n_in();
int {{ con_h_e.name }}_h_e_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ con_h_e.name }}_PHI_E_CONSTRAINT
