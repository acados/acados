#ifndef {{ con_phi_e.name }}_PHI_E_CONSTRAINT
#define {{ con_phi_e.name }}_PHI_E_CONSTRAINT

#ifdef __cplusplus
extern "C" {
#endif

{% if dims.nphi_e > 0 %}
int {{ con_phi_e.name }}_phi_e_constraint(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int {{ con_phi_e.name }}_phi_e_constraint_work(int *, int *, int *, int *);
const int *{{ con_phi_e.name }}_phi_e_constraint_sparsity_in(int);
const int *{{ con_phi_e.name }}_phi_e_constraint_sparsity_out(int);
int {{ con_phi_e.name }}_phi_e_constraint_n_in();
int {{ con_phi_e.name }}_phi_e_constraint_n_out();
{% endif %}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // {{ con_phi_e.name }}_PHI_E_CONSTRAINT
