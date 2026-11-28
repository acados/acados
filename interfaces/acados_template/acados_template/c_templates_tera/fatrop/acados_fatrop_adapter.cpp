/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

// fatrop backend for acados OCPs: fatrop's Nlp<OcpType> on top of an acados ocp_nlp_solver.
//
// All evaluations go through acados' NLP modules (dynamics, cost, constraints), the same code
// paths as acados' SQP: ocp_nlp_approximate_qp_matrices for derivatives, the modules'
// compute_fun for function values. This file only maps acados' QP approximation (qp_in: BAbt,
// RSQrq, DCt, idxb) and function values to fatrop's stage-wise interface.
//
// Constraint rows of a stage are acados' [bounds (idxb); g; h]. Rows acados marks as equalities
// (idxe, e.g. the initial state) are fatrop equalities, rows without finite bound are dropped, all
// others are inequalities. fatrop has one signed multiplier m per row, acados two nonnegative
// ones: lam_lower = max(-m, 0), lam_upper = max(m, 0).
//
// fatrop must be built against acados' BLASFEO: matrices are passed between both directly.

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fatrop/fatrop.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

extern "C"
{
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_constraints_bgh.h"
#include "acados_c/ocp_nlp_interface.h"
}

namespace
{
using fatrop::Hessian;
using fatrop::Index;
using fatrop::Jacobian;
using fatrop::NlpDims;
using fatrop::OcpType;
using fatrop::ProblemDims;
using fatrop::Scalar;
using fatrop::VecRealView;
typedef fatrop::ProblemInfo<OcpType> OcpInfo;

/// Constraint rows of one stage and their fatrop column.
struct StageRows
{
    int nu = 0, nx = 0, nb = 0, ng = 0, nh = 0;
    const int* idxb = nullptr;            // bound row i acts on [u; x][idxb[i]]
    const blasfeo_dvec* d = nullptr;      // acados bounds [lb; ub]
    const blasfeo_dvec* dmask = nullptr;  // 0 for an infinite bound
    std::vector<char> is_eq;              // equality rows of acados (idxe)
    // row r -> fatrop equality column (>= 0), inequality column (-target-1) or FREE
    static constexpr int FREE = INT_MIN;
    std::vector<int> target;
    int neq = 0, nineq = 0;

    int ni() const { return nb + ng + nh; }
    int nv() const { return nu + nx; }
    double lb(int r) const { return BLASFEO_DVECEL(dmask, r) ? BLASFEO_DVECEL(d, r) : -INFINITY; }
    double ub(int r) const
    {
        return BLASFEO_DVECEL(dmask, ni() + r) ? BLASFEO_DVECEL(d, ni() + r) : INFINITY;
    }
    /// Row type: 0 (acados equality), FREE (no finite bound) or -1 (inequality).
    int type(int r) const
    {
        if (is_eq[r]) return 0;
        return std::isinf(lb(r)) && std::isinf(ub(r)) ? FREE : -1;
    }
    /// Signed multiplier of row r from fatrop's equality / inequality multipliers.
    double mult(int r, const double* eq, const double* ineq) const
    {
        const int t = target[r];
        return t >= 0 ? eq[t] : t == FREE ? 0. : ineq[-t - 1];
    }
};

class AcadosFatropNlp : public fatrop::Nlp<OcpType>
{
 public:
    AcadosFatropNlp(ocp_nlp_solver* solver, ocp_nlp_in* in, ocp_nlp_out* out)
        : config_(solver->config),
          dims_(solver->dims),
          in_(in),
          out_(out),
          N_(dims_->N),
          rows_(make_rows(dims_, in_)),
          ocp_dims_(make_dims(rows_)),
          nlp_dims_(make_nlp_dims(ocp_dims_))
    {
        config_->get(config_, dims_, solver->mem, "nlp_mem", &mem_);
        config_->opts_get(config_, solver->opts, "nlp_opts", &opts_);
        config_->work_get(config_, dims_, solver->work, "nlp_work", &work_);
    }

    /// Prepare a solve: initialize the acados modules (as acados' SQP does) and take the initial
    /// guess from nlp_out. Returns false if an inequality row changed between bounded and free.
    bool prepare()
    {
        for (const StageRows& s : rows_)
            for (int r = 0; r < s.ni(); ++r)
                if (!s.is_eq[r] &&
                    (s.type(r) == StageRows::FREE) != (s.target[r] == StageRows::FREE))
                    return false;
        ocp_nlp_initialize_submodules(config_, dims_, in_, out_, opts_, mem_, work_);
        init_ux_.resize(N_ + 1);
        for (int k = 0; k <= N_; ++k)
        {
            init_ux_[k].resize(rows_[k].nv());
            blasfeo_unpack_dvec(rows_[k].nv(), out_->ux + k, 0, init_ux_[k].data(), 1);
        }
        lin_valid_ = false;
        return true;
    }

    /// Write primal solution and multipliers of the iterate to nlp_out.
    void store_solution(const OcpInfo& info, const VecRealView& primal_x, const VecRealView& lam)
    {
        set_primal(info, primal_x.data());
        set_duals(info, lam.data(), 1.);
    }

    const NlpDims& nlp_dims() const override { return nlp_dims_; }
    const ProblemDims<OcpType>& problem_dims() const override { return ocp_dims_; }

    Index eval_derivatives(const OcpInfo& info, const Scalar objective_scale,
                           const VecRealView& primal_x, const VecRealView& primal_s,
                           const VecRealView& lam, VecRealView& grad_x, VecRealView& grad_s,
                           Jacobian<OcpType>& jac, Hessian<OcpType>& hess) override
    {
        // gradient and Jacobians do not depend on the multipliers: read them from the
        // linearization of the Hessian
        lag_hess(info, objective_scale, primal_x, lam, hess);
        constr_jac(jac);
        objective_gradient(info, objective_scale, grad_x, grad_s);
        return 0;
    }

    Index eval_lag_hess(const OcpInfo& info, const Scalar objective_scale,
                        const VecRealView& primal_x, const VecRealView& primal_s,
                        const VecRealView& lam, Hessian<OcpType>& hess) override
    {
        lag_hess(info, objective_scale, primal_x, lam, hess);
        return 0;
    }

    Index eval_constr_jac(const OcpInfo& info, const VecRealView& primal_x,
                          const VecRealView& primal_s, Jacobian<OcpType>& jac) override
    {
        linearize(info, primal_x, nullptr, 0.);
        constr_jac(jac);
        return 0;
    }

    Index eval_constraint_violation(const OcpInfo& info, const VecRealView& primal_x,
                                    const VecRealView& primal_s, VecRealView& res) override
    {
        set_primal(info, primal_x.data());
        Scalar* rp = res.data();
        for (int k = 0; k < N_; ++k)
        {
            // phi(x, u) - x_next
            config_->dynamics[k]->compute_fun(config_->dynamics[k], dims_->dynamics[k],
                                              in_->dynamics[k], opts_->dynamics[k],
                                              mem_->dynamics[k], work_->dynamics[k]);
            blasfeo_unpack_dvec(dims_->nx[k + 1],
                                config_->dynamics[k]->memory_get_fun_ptr(mem_->dynamics[k]), 0,
                                rp + info.offsets_g_eq_dyn[k], 1);
        }
        ocp_nlp_eval_constraints_common(config_, dims_, in_, out_, opts_, mem_, work_);
        for (int k = 0; k <= N_; ++k)
        {
            const StageRows& s = rows_[k];
            // fun = [lb - g; g - ub], masked
            const blasfeo_dvec* fun = mem_->ineq_fun + k;
            double* geq = rp + info.offsets_g_eq_path[k];
            double* gineq = rp + info.offsets_g_eq_slack[k];
            for (int r = 0; r < s.ni(); ++r)
            {
                const int t = s.target[r];
                const double l = s.lb(r);
                const double g = std::isinf(l) ? BLASFEO_DVECEL(fun, s.ni() + r) + s.ub(r) :
                                                 l - BLASFEO_DVECEL(fun, r);
                if (t >= 0)
                    geq[t] = g - l;
                else if (t != StageRows::FREE)
                    gineq[-t - 1] = g;
            }
        }
        res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) =
            res.block(info.number_of_g_eq_slack, info.offset_g_eq_slack) -
            primal_s.block(info.number_of_g_eq_slack, 0);
        return 0;
    }

    Index eval_objective_gradient(const OcpInfo& info, const Scalar objective_scale,
                                  const VecRealView& primal_x, const VecRealView& primal_s,
                                  VecRealView& grad_x, VecRealView& grad_s) override
    {
        linearize(info, primal_x, nullptr, 0.);
        objective_gradient(info, objective_scale, grad_x, grad_s);
        return 0;
    }

    Index eval_objective(const OcpInfo& info, const Scalar objective_scale,
                         const VecRealView& primal_x, const VecRealView& primal_s,
                         Scalar& res) override
    {
        set_primal(info, primal_x.data());
        ocp_nlp_cost_compute(config_, dims_, in_, out_, opts_, mem_, work_);
        res = objective_scale * mem_->cost_value;
        return 0;
    }

    Index get_bounds(const OcpInfo& info, VecRealView& lower_bounds,
                     VecRealView& upper_bounds) override
    {
        for (int k = 0; k <= N_; ++k)
        {
            const StageRows& s = rows_[k];
            for (int r = 0; r < s.ni(); ++r)
            {
                const int t = s.target[r];
                if (t >= 0 || t == StageRows::FREE) continue;
                lower_bounds.data()[info.offsets_slack[k] - t - 1] = s.lb(r);
                upper_bounds.data()[info.offsets_slack[k] - t - 1] = s.ub(r);
            }
        }
        return 0;
    }

    Index get_initial_primal(const OcpInfo& info, VecRealView& primal_x) override
    {
        Scalar* xp = primal_x.data();
        for (int k = 0; k <= N_; ++k)
        {
            const std::vector<double>& ux = init_ux_[k];
            std::copy(ux.begin(), ux.begin() + rows_[k].nu, xp + info.offsets_primal_u[k]);
            std::copy(ux.begin() + rows_[k].nu, ux.end(), xp + info.offsets_primal_x[k]);
        }
        return 0;
    }

    void get_primal_damping(const OcpInfo& info, VecRealView& damping) override { damping = 0; }

    void apply_jacobian_s_transpose(const OcpInfo& info, const VecRealView& multipliers,
                                    const Scalar alpha, const VecRealView& y,
                                    VecRealView& out) override
    {
        out = alpha * y;
        out.block(info.number_of_slack_variables, 0) =
            out.block(info.number_of_slack_variables, 0) -
            multipliers.block(info.number_of_slack_variables, info.offset_g_eq_slack);
    }

 private:
    /// Constraint rows of every stage; equality / inequality / free from the current bounds.
    static std::vector<StageRows> make_rows(ocp_nlp_dims* dims, ocp_nlp_in* in)
    {
        std::vector<StageRows> rows(dims->N + 1);
        for (int k = 0; k <= dims->N; ++k)
        {
            auto* cdims = static_cast<ocp_nlp_constraints_bgh_dims*>(dims->constraints[k]);
            auto* model = static_cast<ocp_nlp_constraints_bgh_model*>(in->constraints[k]);
            StageRows& s = rows[k];
            s.nu = dims->nu[k];
            s.nx = dims->nx[k];
            s.nb = cdims->nb;
            s.ng = cdims->ng;
            s.nh = cdims->nh;
            s.idxb = model->idxb;
            s.d = &model->d;
            s.dmask = model->dmask;
            s.is_eq.assign(s.ni(), 0);
            for (int i = 0; i < cdims->nbue + cdims->nbxe + cdims->nge + cdims->nhe; ++i)
                s.is_eq[model->idxe[i]] = 1;
            s.target.resize(s.ni());
            for (int r = 0; r < s.ni(); ++r)
            {
                const int t = s.type(r);
                s.target[r] = t == StageRows::FREE ? t : t == 0 ? s.neq++ : -(++s.nineq);
            }
        }
        return rows;
    }
    static ProblemDims<OcpType> make_dims(const std::vector<StageRows>& rows)
    {
        const Index K = rows.size();
        std::vector<Index> nu(K), nx(K), neq(K), nineq(K);
        for (Index k = 0; k < K; ++k)
        {
            nu[k] = rows[k].nu;
            nx[k] = rows[k].nx;
            neq[k] = rows[k].neq;
            nineq[k] = rows[k].nineq;
        }
        return ProblemDims<OcpType>(K, nu, nx, neq, nineq);
    }
    static NlpDims make_nlp_dims(const ProblemDims<OcpType>& d)
    {
        Index nvar = 0, neq = 0, nineq = 0;
        for (Index k = 0; k < d.K; ++k)
        {
            nvar += d.number_of_controls[k] + d.number_of_states[k];
            neq += d.number_of_eq_constraints[k] + d.number_of_ineq_constraints[k];
            if (k != d.K - 1) neq += d.number_of_states[k + 1];
            nineq += d.number_of_ineq_constraints[k];
        }
        return NlpDims(nvar, neq, nineq);
    }

    /// nlp_out: u, x = primal iterate (written directly, as acados' own solvers do).
    void set_primal(const OcpInfo& info, const double* xp)
    {
        double* p = const_cast<double*>(xp);
        for (int k = 0; k <= N_; ++k)
        {
            blasfeo_pack_dvec(rows_[k].nu, p + info.offsets_primal_u[k], 1, out_->ux + k, 0);
            blasfeo_pack_dvec(rows_[k].nx, p + info.offsets_primal_x[k], 1, out_->ux + k,
                              rows_[k].nu);
        }
    }

    /// nlp_out: pi, lam = scale * multipliers of fatrop's iterate.
    void set_duals(const OcpInfo& info, const double* lam, double scale)
    {
        for (int k = 0; k <= N_; ++k)
        {
            const StageRows& s = rows_[k];
            if (k < N_)
                for (int i = 0; i < dims_->nx[k + 1]; ++i)
                    BLASFEO_DVECEL(out_->pi + k, i) = scale * lam[info.offsets_g_eq_dyn[k] + i];
            const double* eq = lam + info.offsets_g_eq_path[k];
            const double* ineq = lam + info.offsets_g_eq_slack[k];
            for (int r = 0; r < s.ni(); ++r)
            {
                const double m = scale * s.mult(r, eq, ineq);
                BLASFEO_DVECEL(out_->lam + k, r) = std::max(-m, 0.);
                BLASFEO_DVECEL(out_->lam + k, s.ni() + r) = std::max(m, 0.);
            }
            // no multipliers on infinite bounds, as ocp_nlp_out_set
            blasfeo_dvecmul(2 * s.ni(), in_->dmask + k, 0, out_->lam + k, 0, out_->lam + k, 0);
        }
    }

    /// Hessian of s * f + lam^T c. acados scales only the cost: s * (f + (lam / s)^T c) for s > 0,
    /// otherwise (f + lam^T c) - f.
    void lag_hess(const OcpInfo& info, double s, const VecRealView& primal_x,
                  const VecRealView& lam, Hessian<OcpType>& hess)
    {
        linearize(info, primal_x, lam.data(), s > 0. ? 1. / s : 1.);
        for (int k = 0; k <= N_; ++k)
            blasfeo_dgecpsc(rows_[k].nv(), rows_[k].nv(), s > 0. ? s : 1., mem_->qp_in->RSQrq + k,
                            0, 0, &hess.RSQrqt[k].mat(), 0, 0);
        if (!(s > 0.))
        {
            linearize(info, primal_x, lam.data(), 0.);
            for (int k = 0; k <= N_; ++k)
                blasfeo_dgead(rows_[k].nv(), rows_[k].nv(), -1., mem_->qp_in->RSQrq + k, 0, 0,
                              &hess.RSQrqt[k].mat(), 0, 0);
        }
        // acados fills the lower triangle, fatrop reads the full matrix
        for (int k = 0; k <= N_; ++k)
            blasfeo_dtrtr_l(rows_[k].nv(), &hess.RSQrqt[k].mat(), 0, 0, &hess.RSQrqt[k].mat(), 0,
                            0);
    }

    /// Constraint Jacobians from the current linearization.
    void constr_jac(Jacobian<OcpType>& jac)
    {
        for (int k = 0; k <= N_; ++k)
        {
            const StageRows& s = rows_[k];
            if (k < N_)
                blasfeo_dgecp(s.nv(), dims_->nx[k + 1], mem_->qp_in->BAbt + k, 0, 0,
                              &jac.BAbt[k].mat(), 0, 0);
            MAT* Geq = &jac.Gg_eqt[k].mat();
            MAT* Gineq = &jac.Gg_ineqt[k].mat();
            if (s.neq) blasfeo_dgese(s.nv(), s.neq, 0., Geq, 0, 0);
            if (s.nineq) blasfeo_dgese(s.nv(), s.nineq, 0., Gineq, 0, 0);
            auto set = [&](int r, int var, double v)
            {
                const int t = s.target[r];
                if (t >= 0)
                    BLASFEO_DMATEL(Geq, var, t) = v;
                else if (t != StageRows::FREE)
                    BLASFEO_DMATEL(Gineq, var, -t - 1) = v;
            };
            for (int i = 0; i < s.nb; ++i) set(i, s.idxb[i], 1.);
            const blasfeo_dmat* DCt = mem_->qp_in->DCt + k;
            for (int j = 0; j < s.ng + s.nh; ++j)
                for (int var = 0; var < s.nv(); ++var)
                    set(s.nb + j, var, BLASFEO_DMATEL(DCt, var, j));
        }
    }

    /// Scaled objective gradient from the current linearization.
    void objective_gradient(const OcpInfo& info, double objective_scale, VecRealView& grad_x,
                            VecRealView& grad_s)
    {
        grad_s.block(info.number_of_g_eq_slack, 0) = 0;
        for (int k = 0; k <= N_; ++k)
            blasfeo_dveccpsc(rows_[k].nv(), objective_scale, mem_->cost_grad + k, 0, &grad_x.vec(),
                             grad_x.ai() + info.offsets_tangent_u[k]);
    }

    /// Evaluate acados' QP approximation at primal_x (and, if lam is given, with the multipliers
    /// lam * lam_scale). Cached: gradient, Jacobians and Hessian of an iterate share it.
    void linearize(const OcpInfo& info, const VecRealView& primal_x, const double* lam,
                   double lam_scale)
    {
        const double* xp = primal_x.data();
        const Index nx_tot = nlp_dims_.number_of_variables;
        const Index nlam = nlp_dims_.number_of_eq_constraints;
        const bool same_x = lin_valid_ && std::equal(xp, xp + nx_tot, lin_x_.begin());
        if (same_x && (!lam || (lin_has_duals_ && lam_scale == lin_lam_scale_ &&
                                std::equal(lam, lam + nlam, lin_lam_.begin()))))
            return;
        set_primal(info, xp);
        if (lam) set_duals(info, lam, lam_scale);
        ocp_nlp_approximate_qp_matrices(config_, dims_, in_, out_, opts_, mem_, work_);
        lin_valid_ = true;
        lin_x_.assign(xp, xp + nx_tot);
        lin_has_duals_ = lam != nullptr;
        if (lam)
        {
            lin_lam_.assign(lam, lam + nlam);
            lin_lam_scale_ = lam_scale;
        }
    }

    ocp_nlp_config* config_;
    ocp_nlp_dims* dims_;
    ocp_nlp_in* in_;
    ocp_nlp_out* out_;
    ocp_nlp_memory* mem_ = nullptr;
    ocp_nlp_opts* opts_ = nullptr;
    ocp_nlp_workspace* work_ = nullptr;
    int N_;
    std::vector<StageRows> rows_;
    ProblemDims<OcpType> ocp_dims_;
    NlpDims nlp_dims_;
    std::vector<std::vector<double>> init_ux_;
    // linearization cache
    bool lin_valid_ = false, lin_has_duals_ = false;
    double lin_lam_scale_ = 0.;
    std::vector<double> lin_x_, lin_lam_;
};

struct Solver
{
    std::shared_ptr<AcadosFatropNlp> nlp;
    fatrop::OptionRegistry options;
    std::unique_ptr<fatrop::IpAlgBuilder<OcpType>> builder;
    std::shared_ptr<fatrop::IpAlgorithm<OcpType>> alg;
    std::shared_ptr<fatrop::IpData<OcpType>> data;
    int status = -1;
    double time_total = 0.;
};

/// Runs f, catching C++ exceptions at the C boundary; returns 0 on success.
template <typename F>
int guarded(const char* fn, F&& f)
{
    try
    {
        f();
        return 0;
    }
    catch (const std::exception& e)
    {
        printf("\nerror: acados_fatrop_%s: %s\n", fn, e.what());
        return 1;
    }
}
}  // namespace

#ifdef _WIN32
#define ACADOS_FATROP_API __declspec(dllexport)
#else
#define ACADOS_FATROP_API __attribute__((visibility("default")))
#endif

extern "C"
{
    /// fatrop solver on top of an acados solver (ocp_nlp_solver, nlp_in, nlp_out). Inequality rows
    /// without finite bound at this point are dropped.
    ACADOS_FATROP_API void* acados_fatrop_create(void* nlp_solver, void* nlp_in, void* nlp_out)
    {
        try
        {
            std::unique_ptr<Solver> s(new Solver());
            s->nlp = std::make_shared<AcadosFatropNlp>(static_cast<ocp_nlp_solver*>(nlp_solver),
                                                       static_cast<ocp_nlp_in*>(nlp_in),
                                                       static_cast<ocp_nlp_out*>(nlp_out));
            s->builder.reset(new fatrop::IpAlgBuilder<OcpType>(s->nlp));
            s->alg = s->builder->with_options_registry(&s->options).build();
            s->data = s->builder->get_ipdata();
            return s.release();
        }
        catch (const std::exception& e)
        {
            printf("\nerror: acados_fatrop_create: %s\n", e.what());
            return nullptr;
        }
    }

    ACADOS_FATROP_API void acados_fatrop_free(void* h) { delete static_cast<Solver*>(h); }

    /// Set a fatrop option to value, or to str if not null, converted to the option's type; returns
    /// 2 for an unknown name, 3 for a value of the wrong type and 1 on other errors.
    ACADOS_FATROP_API int acados_fatrop_option_set(void* h, const char* name, double value,
                                                   const char* str)
    {
        using fatrop::OptionTypeId;
        fatrop::OptionRegistry& opts = static_cast<Solver*>(h)->options;
        const OptionTypeId type = opts.get_option_type(name);
        const bool valid = type == OptionTypeId::String ?
                               str != nullptr :
                               !str &&
                                   (type == OptionTypeId::Scalar || value == std::floor(value)) &&
                                   (type != OptionTypeId::Bool || value == 0. || value == 1.);
        if (type == OptionTypeId::Unknown || !valid)
        {
            printf("\nerror: acados_fatrop_option_set: %s option %s\n",
                   type == OptionTypeId::Unknown ? "unknown" : "wrong type for", name);
            return type == OptionTypeId::Unknown ? 2 : 3;
        }
        return guarded("option_set",
                       [&]
                       {
                           if (type == OptionTypeId::String)
                               opts.set_option(name, std::string(str));
                           else if (type == OptionTypeId::Scalar)
                               opts.set_option(name, fatrop::Scalar(value));
                           else if (type == OptionTypeId::Index)
                               opts.set_option(name, fatrop::Index(value));
                           else
                               opts.set_option(name, value != 0.);
                       });
    }

    /// Solve from the iterate in nlp_out and write the result back to it; returns fatrop's
    /// IpSolverReturnFlag, -1 on an exception and -2 if an inequality row changed between bounded
    /// and free (no finite bound) since _create.
    ACADOS_FATROP_API int acados_fatrop_solve(void* h)
    {
        Solver* s = static_cast<Solver*>(h);
        s->status = -1;
        try
        {
            if (!s->nlp->prepare()) return s->status = -2;
            s->status = static_cast<int>(s->alg->optimize());
            const auto& it = s->data->current_iterate();
            s->nlp->store_solution(s->data->info(), it.primal_x(), it.dual_eq());
        }
        catch (const std::exception& e)
        {
            printf("\nerror: acados_fatrop_solve: %s\n", e.what());
            return s->status;
        }
        s->time_total = s->data->timing_statistics().full_algorithm.elapsed();
        return s->status;
    }

    /// Stats: nlp_iter, time_tot [s], status (fatrop's return flag).
    ACADOS_FATROP_API int acados_fatrop_get_stat(void* h, const char* field, double* out)
    {
        return guarded("get_stat",
                       [&]
                       {
                           Solver* s = static_cast<Solver*>(h);
                           if (!strcmp(field, "nlp_iter"))
                               *out = s->data->iteration_number();
                           else if (!strcmp(field, "time_tot"))
                               *out = s->time_total;
                           else if (!strcmp(field, "status"))
                               *out = s->status;
                           else
                               throw std::invalid_argument(std::string("unknown field ") + field);
                       });
    }
}
