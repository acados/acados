
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"

#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"

#include "blasfeo/include/blasfeo_d_aux.h"

namespace acados
{
using std::vector;

ocp_nlp_solution::ocp_nlp_solution(std::shared_ptr<ocp_nlp_out> solution,
                                   std::shared_ptr<ocp_nlp_dims> dims, int status)
    : N(dims->N), nlp_out_(nullptr), dims_(nullptr), status_(status)
{
    if (solution == nullptr || dims == nullptr)
        throw std::invalid_argument("Null pointer passed to constructor");

    nlp_out_ = solution;
    dims_ = dims;
}

ocp_nlp_solution::ocp_nlp_solution(const ocp_nlp_solution &other)
    : N(other.N), nlp_out_(nullptr), dims_(nullptr)
{
    nlp_out_ = other.nlp_out_;
    dims_ = other.dims_;
    status_ = other.status_;
}

vector<vector<double>> ocp_nlp_solution::states()
{
    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage)
    {
        int nx = dims_->nx[stage];
        vector<double> tmp(nx);
        blasfeo_unpack_dvec(nx, &nlp_out_->ux[stage], dims_->nu[stage], tmp.data());
        result.push_back(tmp);
    }
    return result;
}

vector<vector<double>> ocp_nlp_solution::controls()
{
    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage)
    {
        int nu = dims_->nu[stage];
        vector<double> tmp(nu);
        blasfeo_unpack_dvec(nu, &nlp_out_->ux[stage], 0, tmp.data());
        result.push_back(tmp);
    }
    return result;
}

vector<vector<double>> ocp_nlp_solution::lag_mul_bounds()
{
    ocp_nlp_constraints_bgh_dims **constraint_dims = (ocp_nlp_constraints_bgh_dims **)
                                                     dims_->constraints;

    vector<vector<double>> result;
    for (int stage = 0; stage <= N; ++stage)
    {
        int nbx = constraint_dims[stage]->nbx;
        int nbu = constraint_dims[stage]->nbu;

        // extract lam_x
        vector<double> tmp_x(nbx);
        blasfeo_unpack_dvec(nbx, &nlp_out_->pi[stage], nbu, tmp_x.data());

        // extract lam_u
        vector<double> tmp_u(nbu);
        blasfeo_unpack_dvec(nbu, &nlp_out_->pi[stage], 0, tmp_u.data());

        // concatenate
        tmp_x.insert(tmp_x.end(), tmp_u.begin(), tmp_u.end());
        result.push_back(tmp_x);
    }
    return result;
}

vector<vector<double>> ocp_nlp_solution::lag_mul_dynamics()
{
    vector<vector<double>> result;
    for (int stage = 0; stage < N; ++stage)
    {
        int nx1 = dims_->nx[stage + 1];
        vector<double> tmp(nx1);
        blasfeo_unpack_dvec(nx1, &nlp_out_->lam[stage], 0, tmp.data());
        result.push_back(tmp);
    }
    return result;
}

vector<vector<double>> ocp_nlp_solution::lag_mul_constraints()
{
    vector<vector<double>> result;
    return result;
}

ocp_nlp_info ocp_nlp_solution::info() {
    return {
        nlp_out_->sqp_iter,
        status_,
        nlp_out_->inf_norm_res,
        nlp_out_->total_time,
        nlp_out_->qp_iter
    };
}

}  // namespace acados
