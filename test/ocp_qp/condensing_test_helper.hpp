#include "acados/ocp_qp/condensing.h"

void calculate_num_state_bounds(const ocp_qp_in *in, condensing_workspace *work);

int_t get_num_condensed_vars(const ocp_qp_in *in);

int_t get_num_constraints(const ocp_qp_in *in, condensing_workspace *work);

void fill_in_condensing_structs(const ocp_qp_in * const qp_in, condensing_in *in,
    condensing_out *out, condensing_workspace *work);
