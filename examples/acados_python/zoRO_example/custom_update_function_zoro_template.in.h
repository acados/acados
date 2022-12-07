#include "acados_solver_{{ model.name }}.h"

#define {{ model.name | upper }}_ZO_NLBXT  {{zoro_stuff.nlbx_t}}
#define {{ model.name | upper }}_ZO_NUBXT  {{zoro_stuff.nubx_t}}
#define {{ model.name | upper }}_ZO_NLBUT  {{zoro_stuff.nlbu_t}}
#define {{ model.name | upper }}_ZO_NUBUT  {{zoro_stuff.nubu_t}}
#define {{ model.name | upper }}_ZO_NLGT   {{zoro_stuff.nlg_t}}
#define {{ model.name | upper }}_ZO_NUGT   {{zoro_stuff.nug_t}}
#define {{ model.name | upper }}_ZO_NLHT   {{zoro_stuff.nlh_t}}
#define {{ model.name | upper }}_ZO_NUHT   {{zoro_stuff.nuh_t}}
#define {{ model.name | upper }}_ZO_NLBXNT  {{zoro_stuff.nlbx_e_t}}
#define {{ model.name | upper }}_ZO_NUBXNT  {{zoro_stuff.nubx_e_t}}
#define {{ model.name | upper }}_ZO_NLGNT   {{zoro_stuff.nlg_e_t}}
#define {{ model.name | upper }}_ZO_NUGNT   {{zoro_stuff.nug_e_t}}
#define {{ model.name | upper }}_ZO_NLHNT   {{zoro_stuff.nlh_e_t}}
#define {{ model.name | upper }}_ZO_NUHNT   {{zoro_stuff.nuh_e_t}}

// Called at the end of solver creation.
// This is allowed to allocate memory and store the pointer to it into capsule->custom_update_memory.
int custom_update_init_function({{ model.name }}_solver_capsule* capsule);


// Custom update function that can be called between solver calls
int custom_update_function({{ model.name }}_solver_capsule* capsule, double* data, int data_len);


// Called just before destroying the solver.
// Responsible to free allocated memory, stored at capsule->custom_update_memory.
int custom_update_terminate_function({{ model.name }}_solver_capsule* capsule);