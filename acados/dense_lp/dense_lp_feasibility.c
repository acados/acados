/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stddef.h>

// acados
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/dense_lp/dense_lp_feasibility.h"

// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp_utils.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"



#ifndef SIMPLEX_PROFILING
#define SIMPLEX_PROFILING 0
#endif

#ifndef PARTIAL_PRICING
#define PARTIAL_PRICING 1
#endif

#ifndef SIMPLEX_PRICING_BLOCK_SIZE
#define SIMPLEX_PRICING_BLOCK_SIZE 64
#endif

/* If m * N is below this threshold, the GEMV for full structural pricing is
 * cheap enough that all blocks are priced in one shot (no rotation). */
#ifndef SIMPLEX_FULL_PRICING_MAX_MN
#define SIMPLEX_FULL_PRICING_MAX_MN 400000
#endif

/* Scale equality and generic rows to unit 2-norm. Off by default, since it
 * changes the Phase-I objective and doubles the iterations on random LPs.
 * tol_infeas always applies to the row-normalized violation. */
#ifndef SIMPLEX_ROW_NORMALIZATION
#define SIMPLEX_ROW_NORMALIZATION 0
#endif

/* Maximum number of LU updates of the reduced basis before a fresh
 * factorization is forced (bounds the accumulation of rounding errors). */
#ifndef SIMPLEX_LU_UPDATE_MAX
#define SIMPLEX_LU_UPDATE_MAX 64
#endif

/* Debug: check ||P*T - L*U|| after every factor update (O(k^3)). */
#ifndef SIMPLEX_DEBUG_LU_UPDATE
#define SIMPLEX_DEBUG_LU_UPDATE 0
#endif

/*
 * Dual phase for early infeasibility detection. The solve starts with dual
 * simplex iterations (bound-flipping ratio test) on a perturbed cost that
 * makes the slack basis dual feasible:
 *     c_j = +w_j at lower, -w_j at upper, 0 if free or fixed,
 *     w_j = ||a_j|| (1 + u_j / 2), u_j pseudo-random in [0, 1).
 * An empty ratio test yields a Farkas certificate. With a fraction f of the
 * rows violated at the start, the phase runs for at most
 * max(SIMPLEX_DUAL_MIN_STEPS, SIMPLEX_DUAL_ITER_BUDGET * f^2 * (m + nv))
 * steps, so a nearly feasible start only gets a short probe. Afterwards, or
 * on numerical trouble, primal Phase I continues from the current basis.
 */
#ifndef SIMPLEX_DUAL_DETECTION
#define SIMPLEX_DUAL_DETECTION 1
#endif

/* dual steps for a start violating every row, as a fraction of m + nv
 * (scaled by the squared fraction of rows violated at the start) */
#ifndef SIMPLEX_DUAL_ITER_BUDGET
#define SIMPLEX_DUAL_ITER_BUDGET 0.4
#endif

/* minimum number of dual steps (probe for easy certificates) */
#ifndef SIMPLEX_DUAL_MIN_STEPS
#define SIMPLEX_DUAL_MIN_STEPS 8
#endif

/* relative pivot tolerance on the row entries alpha_j in the dual ratio test */
#ifndef SIMPLEX_DUAL_PIVOT_TOL
#define SIMPLEX_DUAL_PIVOT_TOL 1e-7
#endif

/*
 * Devex pricing (Harris 1973; update of Forrest & Goldfarb 1992) in primal
 * Phase I: candidates are ranked by score_j^2 / w_j. With pivot row
 * alpha_j = (B^{-1} a_j)_r and entering q, the weights are updated as
 *     w_j = max(w_j, (alpha_j / alpha_q)^2 w_q),  w_leaving = max(w_q / alpha_q^2, 1),
 * and reset to 1 at the start, after the dual phase, and above
 * SIMPLEX_DEVEX_RESET. Off by default (Dantzig pricing on the Phase-I slope):
 * with the bound-flipping ratio test, Devex did not reduce the runtime.
 */
#ifndef SIMPLEX_DEVEX
#define SIMPLEX_DEVEX 0
#endif

#ifndef SIMPLEX_DEVEX_RESET
#define SIMPLEX_DEVEX_RESET 1e6
#endif

/* Start of the structural variables. 1: projection of 0 onto [l, u], i.e.,
 * superbasic if l < 0 < u. 0: at a bound, the lower one if finite. */
#ifndef SIMPLEX_SUPERBASIC_START
#define SIMPLEX_SUPERBASIC_START 1
#endif

#ifndef SIMPLEX_RHS_REFRESH_INTERVAL
#define SIMPLEX_RHS_REFRESH_INTERVAL 64
#endif

#ifndef SIMPLEX_CERTIFICATE_TOL
#define SIMPLEX_CERTIFICATE_TOL 1e-7
#endif

/*
 * Early infeasibility detection. Any lambda yields the weak-duality bound
 *     phi* >= lambda^T rhs - sum_j sup_{l_j <= x_j <= u_j} (Aphi^T lambda)_j x_j,
 * so a positive right-hand side proves infeasibility (Farkas). Screened are
 * 1) the Phase-I duals (B^T lambda = c_B), every SIMPLEX_FARKAS_DUAL_INTERVAL
 *    iterations (O(N) with full pricing, which already computes Aphi^T lambda);
 * 2) the row lambda = -+B^{-T} e_r of the most violated basic variable, every
 *    SIMPLEX_FARKAS_ROW_INTERVAL iterations (one dual solve and one GEMV).
 * A passed screen triggers the full rounding-aware check; a failed check
 * postpones the next one by SIMPLEX_FARKAS_BACKOFF iterations. An interval
 * of 0 disables the test.
 */
#ifndef SIMPLEX_FARKAS_DUAL_INTERVAL
#define SIMPLEX_FARKAS_DUAL_INTERVAL 1
#endif

#ifndef SIMPLEX_FARKAS_DUAL_INTERVAL_PARTIAL
#define SIMPLEX_FARKAS_DUAL_INTERVAL_PARTIAL 10
#endif

#ifndef SIMPLEX_FARKAS_ROW_INTERVAL
#define SIMPLEX_FARKAS_ROW_INTERVAL 10
#endif

#ifndef SIMPLEX_FARKAS_BACKOFF
#define SIMPLEX_FARKAS_BACKOFF 10
#endif

struct d_simplex_threshold
{
    double t;
    double alpha;
    int row;
    int side;
};

enum d_simplex_var_status
{
    D_SIMPLEX_BASIC,
    D_SIMPLEX_AT_LOWER,
    D_SIMPLEX_AT_UPPER,
    D_SIMPLEX_FREE_NONBASIC,
    D_SIMPLEX_FIXED_NONBASIC,
    /* nonbasic strictly between its bounds (at least one of them finite);
     * may move in both directions until it reaches a bound or enters */
    D_SIMPLEX_SUPERBASIC
};

enum d_simplex_event_kind
{
    D_SIMPLEX_BASIS_EXCHANGE,
    D_SIMPLEX_BOUND_FLIP,
    D_SIMPLEX_EVENT_FAILURE
};

enum d_simplex_certificate_result
{
    D_SIMPLEX_CERTIFICATE_VALID,
    D_SIMPLEX_CERTIFICATE_NOT_FOUND,
    D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN
};

enum d_simplex_certificate_failure
{
    D_SIMPLEX_CERTIFICATE_FAILURE_NONE,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_INPUT,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_PRODUCT,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_ACCUMULATION,
    D_SIMPLEX_CERTIFICATE_FAILURE_COEFFICIENT_OVERFLOW,
    D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_COEFFICIENT,
    D_SIMPLEX_CERTIFICATE_FAILURE_INFINITE_SUPPORT,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_BOUND,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_SUPPORT,
    D_SIMPLEX_CERTIFICATE_FAILURE_NONPOSITIVE_MARGIN,
    D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_MARGIN,
    D_SIMPLEX_CERTIFICATE_FAILURE_DUAL_RESIDUAL
};

struct d_simplex_certificate_diagnostic
{
    enum d_simplex_certificate_failure failure;
    int index;
    int index_is_variable;
    double value;
    double scale;
    double threshold;
};

struct d_simplex_event
{
    enum d_simplex_event_kind kind;
    int leaving_pos;
    int leaving_side;
    double step;
    double entering_value;
};

// Private internal simplex workspace.
struct d_simplex_workspace
{
    /* original QP dimensions */
    int nv, ne, ng, nb;
    int nv_max, ne_max, ng_max, nb_max, ns_max;

    /* max capacities derived from dims */
    int m_max;
    int nslack_max;
    int N_max;

    /* actual sizes for current solve */
    int m;
    int nslack;
    int nineq_slack;
    int neq_residual;
    int N;

    int initialized;
    int reduced_active;
    int reduced_rep_valid;
    int reduced_lu_valid;
    int nsingleton;
    int nreduced;
    int reduced_k_min;
    int reduced_k_max;
    int reduced_k_samples;
    double reduced_k_sum;
    int pivot_counts[4];
    int bound_flip_count;
    int zero_step_count;
    int certificate_valid;
    double certificate_margin;
    double certificate_scale;
    double certificate_dual_residual;
    double certificate_dual_residual_scale;
    int certificate_check_count;
    int farkas_dual_screen_count;
    int farkas_row_screen_count;
    int farkas_next_full_check_iter;
    int farkas_detection; /* 0: terminal, 1: phase-I duals, 2: single row, 3: dual simplex */
    enum d_simplex_certificate_result certificate_result;
    struct d_simplex_certificate_diagnostic certificate_diagnostic;

    /* basis bookkeeping */
    int *basis_idx;
    int *nonbasis_idx;
    int *ipiv;
    int *row_perm;        /* m_max */
    int *singleton_row;   /* N_max, -1 for non-singleton columns */
    double *singleton_value; /* N_max */
    int *singleton_basis_pos; /* m_max */
    int *singleton_rows;     /* m_max */
    int *reduced_basis_pos;  /* m_max */
    int *reduced_rows;       /* m_max */
    int *row_to_singleton;   /* m_max */
    int *basis_pos_to_reduced; /* m_max */
    int *reduced_var;        /* nv: reduced column -> structural variable */
    int *var_to_reduced;     /* N_max: variable -> reduced column */
    int *row_to_reduced;     /* m_max: physical row -> reduced row */
    int *lambda_nz;          /* m_max: rows with lambda != 0 (if lambda_nnz >= 0) */
    int lambda_nnz;          /* -1: lambda pattern unknown (treat as dense) */
    int *sing_nz;            /* m_max: scratch, singleton positions with lambda != 0 */

    /* LP data / factorization */
    struct blasfeo_dmat Aphi;   /* [m_max x nv_max]: structural columns only, the
                                 * auxiliary columns are singletons (see d_simplex_A) */
    struct blasfeo_dmat B_mat;  /* max-size [m_max x m_max] */
    struct blasfeo_dmat LU_mat;  /* max-size [m_max x m_max] */
    struct blasfeo_dmat G_mat;  /* max-size [m_max x nv] */
    /* Factors of the transposed reduced basis T = H^T: P*T = L*U with L unit
     * lower triangular, U upper triangular, (P*x)[i] = x[row_perm[i]].
     * Rows of T follow reduced_var, columns follow reduced_rows. */
    struct blasfeo_dmat L_red;  /* max-size [nv x nv] */
    struct blasfeo_dmat U_red;  /* max-size [nv x nv] */
    struct blasfeo_dvec rhs;     /* max-size [m_max] */
    struct blasfeo_dvec c;       /* max-size [N_max] */
    struct blasfeo_dvec rc;  /* N_max: c - Aphi^T lambda */
    struct blasfeo_dvec tmp_m_1; /* max-size [m_max] */
    struct blasfeo_dvec tmp_m_2; /* max-size [m_max] */
    struct blasfeo_dvec lu_work;

    double *effective_rhs_data; // m
    double *x_B;        // m
    double *lambda;     // m
    double *d_dir;      // m
    double *cB;         // m
    struct d_simplex_threshold *thresholds; // 2*m_max bound events
    double *w;          // nv
    double *lower;      // N_max
    double *upper;      // N_max
    double *nonbasic_value; // N_max
    double *feas_tol;   // N_max: primal feasibility band of each variable
    double *dual_cost;  // N_max: perturbed cost of the dual simplex phase
    double *dual_d;     // N_max: reduced costs of the dual simplex phase
    double *row_norm;   // m_max: 2-norm of the structural part of each row
    double *devex_w;    // N_max: Devex reference weights
    double *lo_eff;     // N_max: lower bound, -inf if none (ratio test)
    double *up_eff;     // N_max: upper bound, +inf if none (ratio test)
    double *col_norm_sq; // N_max: squared 2-norm of each structural column of Aphi
    double *cert_q_scale; // N_max: scratch, sum_i |Aphi_ij lambda_i| in the certificate check
    unsigned char *has_lower; // N_max
    unsigned char *has_upper; // N_max
    unsigned char *status; // N_max

    int print_level;
    int max_iter;
    double tol;
    int iter;
    double tol_infeas;
    double tol_pivot;
#if SIMPLEX_PROFILING
    double time_rhs_rebuild;
    double time_rhs_update;
    double time_bound_flip;
    double time_recovery_validation;
    double time_factorization;
    double time_basic_solve;
    double time_phase1_objective;
    double time_dual_solve;
    double time_pricing;
    double time_direction_solve;
    double time_ratio_test;
    double time_basis_update;
    double time_reduced_maps;
    double time_reduced_H_build;
    double time_reduced_G_build;
    double time_reduced_factorization_only;
    double time_reduced_incremental_update;
    double time_reduced_triangular_solves;
    double time_reduced_gemv_recovery;
    double time_farkas;
#endif
    int rhs_rebuild_count;
    int rhs_update_count;
    int numerical_refresh_count;
    int rhs_fresh;
    int xB_valid;          /* x_B consistent with the basis and effective rhs */
    int xB_solve_count;
    int events_since_rhs_refresh;
    int numerical_retry_count;
    int reduced_H_scalar_reads;
    int reduced_G_scalar_reads;
    int lu_updates_since_refactor;
    int lu_update_count;
    int lu_update_fail_count;
    int reduced_refactor_count;
    int dual_mode;
    int dual_iter_count;
    int dual_flip_count;
    int dual_fallback_reason; /* 0: none, see d_simplex_dual_fallback_name */
    int dual_detected_infeasible;
    int dual_step_count;      /* dual steps (pivots and flip-only steps) */
    int dual_budget;          /* -1 until set from the start's violated rows */
#if SIMPLEX_PROFILING
    double time_dual_simplex;
    double time_devex;
    int devex_resets;
#endif
};

static const char *d_simplex_certificate_failure_name(
    enum d_simplex_certificate_failure failure)
{
    switch (failure)
    {
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONE:
            return "none";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_INPUT:
            return "non-finite input";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_PRODUCT:
            return "non-finite product";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_ACCUMULATION:
            return "non-finite accumulation";
        case D_SIMPLEX_CERTIFICATE_FAILURE_COEFFICIENT_OVERFLOW:
            return "coefficient overflow";
        case D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_COEFFICIENT:
            return "ambiguous coefficient sign";
        case D_SIMPLEX_CERTIFICATE_FAILURE_INFINITE_SUPPORT:
            return "infinite support";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_BOUND:
            return "non-finite selected bound";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_SUPPORT:
            return "non-finite support accumulation";
        case D_SIMPLEX_CERTIFICATE_FAILURE_NONPOSITIVE_MARGIN:
            return "nonpositive margin";
        case D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_MARGIN:
            return "margin below acceptance threshold";
        case D_SIMPLEX_CERTIFICATE_FAILURE_DUAL_RESIDUAL:
            return "dual residual above tolerance";
    }
    return "unknown";
}

static enum d_simplex_certificate_result d_simplex_certificate_fail(
    struct d_simplex_certificate_diagnostic *diagnostic,
    enum d_simplex_certificate_failure failure,
    int index,
    int index_is_variable,
    double value,
    double scale,
    double threshold,
    enum d_simplex_certificate_result result)
{
    diagnostic->failure = failure;
    diagnostic->index = index;
    diagnostic->index_is_variable = index_is_variable;
    diagnostic->value = (double) value;
    diagnostic->scale = (double) scale;
    diagnostic->threshold = (double) threshold;
    return result;
}

/*
 * rc[off:off+n] = alpha * Aphi[:, off:off+n]^T lambda. If the nonzero
 * pattern of lambda is known and sparse (lambda = B^{-T} c_B is nonzero only
 * on the reduced rows and the rows of violated basic slacks), only those rows
 * of Aphi are visited; otherwise a dense GEMV is used (lambda packed in
 * lu_work by the caller or here).
 */
static void d_simplex_lambda_times_A(struct d_simplex_workspace *ws, double alpha,
                                     int off, int n, int lambda_packed)
{
    int m = ws->m;
    int nnz = ws->lambda_nnz;
    if (n <= 0)
        return;
    if (nnz >= 0 && 4 * nnz < m)
    {
        double *out = ws->rc.pa + off;
        for (int j = 0; j < n; ++j)
            out[j] = 0.0;
        if (nnz == 0)
            return;
        const ptrdiff_t stride = n > 1 ?
            &BLASFEO_DMATEL(&ws->Aphi, 0, 1) - &BLASFEO_DMATEL(&ws->Aphi, 0, 0) : 0;
        for (int q = 0; q < nnz; ++q)
        {
            int i = ws->lambda_nz[q];
            double a = alpha * ws->lambda[i];
            const double *row = &BLASFEO_DMATEL(&ws->Aphi, i, off);
            for (int j = 0; j < n; ++j)
                out[j] += a * row[j * stride];
        }
        return;
    }
    if (!lambda_packed)
        blasfeo_pack_dvec(m, ws->lambda, 1, &ws->lu_work, 0);
    blasfeo_dgemv_t(m, n, alpha, &ws->Aphi, 0, off, &ws->lu_work, 0,
                    0.0, &ws->c, 0, &ws->rc, off);
}

/* Structural entries of q = Aphi^T lambda are in ws->rc, auxiliary ones are
 * -lambda[row]. */
static inline double d_simplex_column_dot_lambda(const struct d_simplex_workspace *ws, int j)
{
    int row = ws->singleton_row[j];
    return row >= 0 ? ws->singleton_value[j] * ws->lambda[row] : BLASFEO_DVECEL(&ws->rc, j);
}

static void d_simplex_compute_dual_product(struct d_simplex_workspace *ws)
{
    d_simplex_lambda_times_A(ws, 1.0, 0, ws->nv, 0);
}

/*
 * Rounding-aware Farkas test of ws->lambda:
 *
 *     margin = lambda^T rhs - sum_j sup_{l_j <= x_j <= u_j} q_j x_j > 0,
 *     q = Aphi^T lambda.
 *
 * Each |q_j| below its cancellation scale tol * sum_i |Aphi_ij lambda_i| (or
 * below an absolute floor tied to the problem scale) is treated as zero, and
 * every support term is charged its worst-case rounding. Non-finite values
 * propagate into margin / scale and are rejected by the final test.
 */
static enum d_simplex_certificate_result d_simplex_check_infeasibility_certificate(
    struct d_simplex_workspace *ws,
    double *margin_out,
    double *scale_out,
    struct d_simplex_certificate_diagnostic *diagnostic)
{
    const double certificate_tolerance = SIMPLEX_CERTIFICATE_TOL;
    const int m = ws->m, nv = ws->nv, N = ws->N;
    const double *lambda = ws->lambda;
    double lambda_rhs = 0.0;
    double lambda_norm1 = 0.0;
    double support = 0.0;
    double scale = 0.0;
    int coefficient_uncertain = 0;

    *margin_out = NAN;
    *scale_out = NAN;
    diagnostic->failure = D_SIMPLEX_CERTIFICATE_FAILURE_NONE;
    diagnostic->index = -1;
    diagnostic->index_is_variable = 0;
    diagnostic->value = NAN;
    diagnostic->scale = NAN;
    diagnostic->threshold = NAN;

    for (int i = 0; i < m; ++i)
    {
        double term = lambda[i] * BLASFEO_DVECEL(&ws->rhs, i);
        lambda_rhs += term;
        lambda_norm1 += fabs(lambda[i]) * ws->row_norm[i];
        scale += fabs(term);
    }

    /* Absolute floor of the numerical-zero test on q_j: with boxed slacks
     * rhs is zero on the generic rows, so sum |lambda_i rhs_i| can vanish;
     * |q_j| <= sum_i |lambda_i| ||row_i|| sets the scale instead. */
    double q_absolute_floor = certificate_tolerance * fmax(scale, lambda_norm1);

    /* q of the structural columns, and q_scale_j = sum_i |Aphi_ij lambda_i|
     * accumulated row-wise over the nonzero lambda_i only. */
    d_simplex_compute_dual_product(ws);
    double *q_scale = ws->cert_q_scale;
    for (int j = 0; j < nv; ++j)
        q_scale[j] = 0.0;
    const ptrdiff_t stride = nv > 1 ?
        &BLASFEO_DMATEL(&ws->Aphi, 0, 1) - &BLASFEO_DMATEL(&ws->Aphi, 0, 0) : 0;
    for (int i = 0; i < m; ++i)
    {
        double abs_lambda = fabs(lambda[i]);
        if (abs_lambda == 0.0)
            continue;
        const double *row = &BLASFEO_DMATEL(&ws->Aphi, i, 0);
        for (int j = 0; j < nv; ++j)
            q_scale[j] += fabs(row[j * stride]) * abs_lambda;
    }

    for (int j = 0; j < N; ++j)
    {
        double q, q_abs_scale;
        int singleton_row = ws->singleton_row[j];
        if (singleton_row >= 0)
        {
            /* auxiliary column: a single term */
            q = ws->singleton_value[j] * lambda[singleton_row];
            q_abs_scale = fabs(q);
        }
        else
        {
            q = BLASFEO_DVECEL(&ws->rc, j);
            q_abs_scale = q_scale[j];
        }
        double q_tolerance = certificate_tolerance * q_abs_scale;

        if (q_abs_scale == 0.0)
        {
            /* every term is exactly zero, so q must be exactly zero */
            if (q != 0.0 && !coefficient_uncertain)
            {
                d_simplex_certificate_fail(
                    diagnostic, D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_COEFFICIENT,
                    j, 1, q, q_abs_scale, q_absolute_floor,
                    D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN);
                coefficient_uncertain = 1;
            }
            continue;
        }
        if (fabs(q) <= q_tolerance || fabs(q) <= q_absolute_floor)
            continue; /* numerically zero */

        int has_lower = ws->has_lower[j];
        int has_upper = ws->has_upper[j];
        if (has_lower && has_upper)
        {
            double contribution = fmax(q * ws->lower[j], q * ws->upper[j]);
            double uncertainty = q_tolerance * fmax(fabs(ws->lower[j]), fabs(ws->upper[j]));
            support += contribution + uncertainty;
            scale += fabs(contribution) + uncertainty;
        }
        else if (has_lower ? q < 0.0 : (has_upper && q > 0.0))
        {
            double bound = has_lower ? ws->lower[j] : ws->upper[j];
            double support_term = q * bound;
            double uncertainty = q_tolerance * fabs(bound);
            support += support_term + uncertainty;
            scale += fabs(support_term) + 2.0 * uncertainty;
        }
        else
        {
            /* q_j pushes towards a missing bound: unbounded support */
            return d_simplex_certificate_fail(
                diagnostic, D_SIMPLEX_CERTIFICATE_FAILURE_INFINITE_SUPPORT,
                j, 1, q, q_abs_scale, q_absolute_floor,
                D_SIMPLEX_CERTIFICATE_NOT_FOUND);
        }
    }

    if (coefficient_uncertain)
        return D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN;

    double margin = lambda_rhs - support;
    *margin_out = margin;
    *scale_out = scale;
    if (!isfinite(margin) || !isfinite(scale))
        return d_simplex_certificate_fail(
            diagnostic, D_SIMPLEX_CERTIFICATE_FAILURE_NONFINITE_ACCUMULATION,
            -1, 0, margin, scale, NAN,
            D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN);
    if (margin <= 0.0)
        return d_simplex_certificate_fail(
            diagnostic, D_SIMPLEX_CERTIFICATE_FAILURE_NONPOSITIVE_MARGIN,
            -1, 0, margin, scale, 0.0, D_SIMPLEX_CERTIFICATE_NOT_FOUND);

    /* accept only margins above a scale-relative allowance for rounding in
     * the dot products and the support accumulation */
    double acceptance_threshold = certificate_tolerance * scale;
    if (margin <= acceptance_threshold)
        return d_simplex_certificate_fail(
            diagnostic, D_SIMPLEX_CERTIFICATE_FAILURE_AMBIGUOUS_MARGIN,
            -1, 0, margin, scale, acceptance_threshold,
            D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN);

    return D_SIMPLEX_CERTIFICATE_VALID;
}

/* ||B^T lambda - c_B||_inf and its scale max_k |q_k| + |c_B,k| (q from
 * d_simplex_compute_dual_product). */
static void d_simplex_compute_dual_residual(
    struct d_simplex_workspace *ws,
    double *residual_out,
    double *scale_out)
{
    double residual_inf = 0.0;
    double scale_inf = 0.0;

    d_simplex_compute_dual_product(ws);
    for (int k = 0; k < ws->m; ++k)
    {
        int basis_var = ws->basis_idx[k];
        double q_var = d_simplex_column_dot_lambda(ws, basis_var);
        residual_inf = fmax(residual_inf, fabs(q_var - ws->cB[k]));
        scale_inf = fmax(scale_inf, fabs(q_var) + fabs(ws->cB[k]));
    }
    *residual_out = residual_inf;
    *scale_out = scale_inf;
}

static enum d_simplex_certificate_result d_simplex_run_certificate_check(
    struct d_simplex_workspace *ws)
{
    double margin = NAN;
    double scale = NAN;
    enum d_simplex_certificate_result result =
        d_simplex_check_infeasibility_certificate(
            ws, &margin, &scale, &ws->certificate_diagnostic);
    ++ws->certificate_check_count;
    ws->certificate_margin = margin;
    ws->certificate_scale = scale;

    if (result == D_SIMPLEX_CERTIFICATE_VALID)
    {
        /* negated test: also rejects a NaN residual */
        double residual_tolerance = SIMPLEX_CERTIFICATE_TOL * ws->certificate_dual_residual_scale;
        if (!(ws->certificate_dual_residual <= residual_tolerance))
        {
            d_simplex_certificate_fail(
                &ws->certificate_diagnostic,
                D_SIMPLEX_CERTIFICATE_FAILURE_DUAL_RESIDUAL,
                -1, 0, ws->certificate_dual_residual,
                ws->certificate_dual_residual_scale,
                residual_tolerance,
                D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN);
            result = D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN;
        }
    }

    ws->certificate_result = result;
    ws->certificate_valid = result == D_SIMPLEX_CERTIFICATE_VALID;
    return result;
}

/*
 * Cheap screen of ws->lambda as a Farkas certificate (margin as in
 * d_simplex_check_infeasibility_certificate, without rounding bookkeeping).
 * q_j = q_sign * rc_j for the structural columns (pricing stores
 * rc = -Aphi^T lambda, so q_sign = -1) and -lambda_row for the auxiliaries.
 * A positive answer only means that the full check is worth running. The
 * numerical-zero scale is sum_i |lambda_i| ||row_i|| >= |q_j|.
 */
static int d_simplex_farkas_screen(const struct d_simplex_workspace *ws, double q_sign)
{
    double margin = 0.0;
    double scale = 0.0;
    double lambda_norm1 = 0.0;
    double q_bound = 0.0; /* sum_i |lambda_i| ||row_i|| >= |q_j| */
    for (int i = 0; i < ws->m; ++i)
    {
        double term = ws->lambda[i] * BLASFEO_DVECEL(&ws->rhs, i);
        lambda_norm1 += fabs(ws->lambda[i]);
        q_bound += fabs(ws->lambda[i]) * ws->row_norm[i];
        margin += term;
        scale += fabs(term);
    }
    if (!(lambda_norm1 > 0.0))
        return 0;

    double q_zero = SIMPLEX_CERTIFICATE_TOL * q_bound;
    for (int j = 0; j < ws->N; ++j)
    {
        int row = ws->singleton_row[j];
        double q = row >= 0 ? ws->singleton_value[j] * ws->lambda[row] :
            q_sign * BLASFEO_DVECEL(&ws->rc, j);
        if (fabs(q) <= q_zero)
            continue;
        double bound;
        if (q > 0.0)
        {
            if (!ws->has_upper[j])
                return 0; /* unbounded support */
            bound = ws->upper[j];
        }
        else
        {
            if (!ws->has_lower[j])
                return 0;
            bound = ws->lower[j];
        }
        double support = q * bound;
        margin -= support;
        scale += fabs(support);
    }
    return isfinite(margin) && margin > SIMPLEX_CERTIFICATE_TOL * scale;
}

/*
 * Full (rounding-aware) certificate check of the current ws->lambda against
 * ws->cB, run after a positive screen. Returns 1 if the certificate is valid.
 * A failed check postpones the next full check by SIMPLEX_FARKAS_BACKOFF
 * iterations. Overwrites ws->rc.
 */
static int d_simplex_try_farkas_certificate(struct d_simplex_workspace *ws,
                                            int iter, int detection_kind)
{
    if (iter < ws->farkas_next_full_check_iter)
        return 0;
    d_simplex_compute_dual_residual(ws, &ws->certificate_dual_residual,
                                    &ws->certificate_dual_residual_scale);
    if (d_simplex_run_certificate_check(ws) == D_SIMPLEX_CERTIFICATE_VALID)
    {
        ws->farkas_detection = detection_kind;
        return 1;
    }
    ws->farkas_next_full_check_iter = iter + SIMPLEX_FARKAS_BACKOFF;
    return 0;
}

static void d_simplex_print_certificate_diagnostic(
    enum d_simplex_certificate_result result,
    const struct d_simplex_certificate_diagnostic *diagnostic)
{
    const char *result_name =
        result == D_SIMPLEX_CERTIFICATE_VALID ? "valid" :
        result == D_SIMPLEX_CERTIFICATE_NOT_FOUND ? "not found" :
        "numerically uncertain";

    printf("\tCertificate check: %s", result_name);
    if (diagnostic->failure == D_SIMPLEX_CERTIFICATE_FAILURE_NONE)
    {
        printf("\n");
        return;
    }
    printf("; failure=%s",
           d_simplex_certificate_failure_name(diagnostic->failure));
    if (diagnostic->index >= 0)
    {
        if (diagnostic->index_is_variable)
            printf("; variable=%d", diagnostic->index);

        else
            printf("; row=%d", diagnostic->index);
    }

    if (diagnostic->index_is_variable)
        printf("; q=%.17e; cancellation_scale=%.17e; sign_tolerance=%.17e\n",
               diagnostic->value, diagnostic->scale, diagnostic->threshold);
    else
        printf("; value=%.17e; scale=%.17e; threshold=%.17e\n",
               diagnostic->value, diagnostic->scale, diagnostic->threshold);
}

static void print_idx_local(const char *name, int *v, int n)
{
    printf("%s = [", name);
    for (int i = 0; i < n; i++)
        printf(" %d", v[i]);
    printf(" ]\n");
}

/* Strict order of breakpoints: step length, then basis position, then side. */
static inline int d_simplex_threshold_less(const struct d_simplex_threshold *a,
                                           const struct d_simplex_threshold *b)
{
    if (a->t != b->t)
        return a->t < b->t;
    if (a->row != b->row)
        return a->row < b->row;
    return a->side < b->side;
}

/* Binary min-heap on the breakpoints (d_simplex_threshold_less). */
static void d_simplex_threshold_sift_down(struct d_simplex_threshold *heap, int n, int i)
{
    struct d_simplex_threshold x = heap[i];
    for (;;)
    {
        int child = 2 * i + 1;
        if (child >= n)
            break;
        if (child + 1 < n && d_simplex_threshold_less(&heap[child + 1], &heap[child]))
            child++;
        if (!d_simplex_threshold_less(&heap[child], &x))
            break;
        heap[i] = heap[child];
        i = child;
    }
    heap[i] = x;
}

static void d_simplex_threshold_heap_pop(struct d_simplex_threshold *heap, int *n,
                                         struct d_simplex_threshold *top)
{
    *top = heap[0];
    heap[0] = heap[--(*n)];
    if (*n > 0)
        d_simplex_threshold_sift_down(heap, *n, 0);
}

/************************************************
* opts
************************************************/

static void dense_lp_feasibility_opts_set_default(dense_lp_feasibility_opts *opts)
{
    opts->print_level = 0;
    opts->max_iter = 2000;
    /* dual (pricing) tolerance on the Phase-I reduced costs */
    opts->tol = 1e-7;
    /* primal feasibility tolerance: max violation of a bound or of a
     * constraint row normalized to unit 2-norm (same measure as DAQP's
     * primal_tol, whose default is 1e-6) */
    opts->tol_infeas = 1e-6;
    opts->tol_pivot = 1e-12;
}

acados_size_t dense_lp_feasibility_opts_calculate_size(void *config_, void *dims_)
{
    (void) config_;
    (void) dims_;

    acados_size_t size = 0;
    size += sizeof(dense_lp_feasibility_opts);
    size += 8;  // align for safety
    make_int_multiple_of(8, &size);

    return size;
}

void *dense_lp_feasibility_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    (void) config_;
    (void) dims_;

    char *c_ptr = (char *) raw_memory;

    align_char_to(8, &c_ptr);

    dense_lp_feasibility_opts *opts = (dense_lp_feasibility_opts *) c_ptr;
    c_ptr += sizeof(dense_lp_feasibility_opts);

    dense_lp_feasibility_opts_set_default(opts);

    assert((char *) raw_memory + dense_lp_feasibility_opts_calculate_size(config_, dims_) >= c_ptr);

    return (void *) opts;
}

void dense_lp_feasibility_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    (void) config_;
    (void) dims_;

    dense_lp_feasibility_opts *opts = (dense_lp_feasibility_opts *) opts_;
    dense_lp_feasibility_opts_set_default(opts);
}

void dense_lp_feasibility_opts_update(void *config_, void *dims_, void *opts_)
{
    (void) config_;
    (void) dims_;
    (void) opts_;

    return;
}

void dense_lp_feasibility_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    (void) config_;

    dense_lp_feasibility_opts *opts = (dense_lp_feasibility_opts *) opts_;

    if (!strcmp(field, "print_level"))
    {
        opts->print_level = *(int *) value;
    }
    else if (!strcmp(field, "max_iter"))
    {
        opts->max_iter = *(int *) value;
    }
    else if (!strcmp(field, "tol"))
    {
        opts->tol = *(double *) value;
    }
    else if (!strcmp(field, "tol_infeas"))
    {
        opts->tol_infeas = *(double *) value;
    }
    else if (!strcmp(field, "tol_pivot"))
    {
        opts->tol_pivot = *(double *) value;
    }
    else
    {
        printf("dense_lp_feasibility_opts_set: got non-supported field %s\n", field);
        exit(1);
    }
}

void dense_lp_feasibility_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    (void) config_;

    dense_lp_feasibility_opts *opts = (dense_lp_feasibility_opts *) opts_;

    if (!strcmp(field, "print_level"))
    {
        *(int *) value = opts->print_level;
    }
    else if (!strcmp(field, "max_iter"))
    {
        *(int *) value = opts->max_iter;
    }
    else if (!strcmp(field, "tol"))
    {
        *(double *) value = opts->tol;
    }
    else if (!strcmp(field, "tol_infeas"))
    {
        *(double *) value = opts->tol_infeas;
    }
    else if (!strcmp(field, "tol_pivot"))
    {
        *(double *) value = opts->tol_pivot;
    }
    else
    {
        printf("dense_lp_feasibility_opts_get: got non-supported field %s\n", field);
        exit(1);
    }
}

/************************************************
* memory
************************************************/

static void dense_lp_feasibility_get_max_sizes(const dense_qp_dims *dims,
                                               int *m_max,
                                               int *nslack_max,
                                               int *N_max)
{
    int nv = dims->nv;
    int ne = dims->ne;
    int ng = dims->ng;

    /* One residual per equality and one boxed slack per generic constraint. */
    *m_max      = ne + ng;
    *nslack_max = *m_max;
    *N_max      = nv + *nslack_max;
}

static acados_size_t dense_lp_feasibility_workspace_calculate_size_internal(void *dims_)
{
    dense_qp_dims *dims = dims_;

    int m_max, nslack_max, N_max;
    dense_lp_feasibility_get_max_sizes(dims, &m_max, &nslack_max, &N_max);

    acados_size_t size = 0;

    // struct itself
    size += sizeof(struct d_simplex_workspace);

    // threshold substructure
    size += (2 * m_max + dims->nv) * sizeof(struct d_simplex_threshold);

    // plain double work arrays
    size += m_max * sizeof(double);   // effective_rhs_data
    size += m_max * sizeof(double);   // x_B
    size += m_max * sizeof(double);   // lambda
    size += m_max * sizeof(double);   // d_dir
    size += m_max * sizeof(double);   // cB
    size += dims->nv * sizeof(double); // w
    size += 3 * N_max * sizeof(double); // lower, upper, nonbasic_value
    size += N_max * sizeof(double); // singleton_value
    size += N_max * sizeof(double); // feas_tol
    size += 2 * N_max * sizeof(double); // dual_cost, dual_d
    size += 2 * N_max * sizeof(double); // col_norm_sq, cert_q_scale
    size += m_max * sizeof(double); // row_norm
    size += 2 * N_max * sizeof(double); // lo_eff, up_eff
    size += N_max * sizeof(double); // devex_w
    size += 3 * N_max * sizeof(unsigned char); // bound flags and status
    size += 8; // align the following integer arrays

    // int arrays
    size += m_max * sizeof(int);   // basis_idx
    size += N_max * sizeof(int);   // nonbasis_idx
    size += m_max * sizeof(int);   // ipiv (or N, but m is enough)
    size += m_max * sizeof(int);  // row_perm
    size += N_max * sizeof(int);  // singleton_row
    size += 6 * m_max * sizeof(int); // reduced basis maps
    size += dims->nv * sizeof(int); // reduced_var
    size += N_max * sizeof(int); // var_to_reduced
    size += m_max * sizeof(int); // row_to_reduced
    size += 2 * m_max * sizeof(int); // lambda_nz, sing_nz

    // BLASFEO descriptors are embedded in the workspace struct. Reserve
    // explicit worst-case padding before the aligned payload blocks.
    size += 64;

    // blasfeo payload memory
    size += blasfeo_memsize_dmat(m_max, dims->nv); // Aphi, structural columns only
    size += blasfeo_memsize_dmat(m_max, m_max); // B_mat
    size += blasfeo_memsize_dmat(m_max, m_max); // LU_mat
    size += blasfeo_memsize_dmat(m_max, dims->nv); // G_mat
    size += 2 * blasfeo_memsize_dmat(dims->nv, dims->nv); // L_red, U_red

    size += blasfeo_memsize_dvec(m_max); // rhs
    size += blasfeo_memsize_dvec(N_max); // c
    size += blasfeo_memsize_dvec(N_max); // rc
    size += blasfeo_memsize_dvec(m_max); // tmp_m_1
    size += blasfeo_memsize_dvec(m_max); // tmp_m_2
    size += blasfeo_memsize_dvec(m_max); // lu_work

    make_int_multiple_of(8, &size);
    return size;
}

static void *dense_lp_feasibility_workspace_assign_internal(void *dims_,
                                                            void *raw_memory,
                                                            struct d_simplex_workspace **ws_out)
{
    dense_qp_dims *dims = dims_;

    int m_max, nslack_max, N_max;
    dense_lp_feasibility_get_max_sizes(dims, &m_max, &nslack_max, &N_max);

    char *c_ptr = (char *) raw_memory;
    align_char_to(8, &c_ptr);

    struct d_simplex_workspace *ws = (struct d_simplex_workspace *) c_ptr;
    c_ptr += sizeof(struct d_simplex_workspace);

    align_char_to(8, &c_ptr);
    ws->thresholds = (struct d_simplex_threshold *) c_ptr;
    /* 2 bound events per row (primal ratio test) or one per nonbasic
     * variable, N - m = nv (dual ratio test) */
    c_ptr += (2 * m_max + dims->nv) * sizeof(struct d_simplex_threshold);

    /* store capacities and original dims */
    ws->nv = dims->nv;
    ws->ne = dims->ne;
    ws->ng = dims->ng;
    ws->nb = dims->nb;
    ws->nv_max = dims->nv;
    ws->ne_max = dims->ne;
    ws->ng_max = dims->ng;
    ws->nb_max = dims->nb;
    ws->ns_max = dims->ns;

    ws->m_max = m_max;
    ws->nslack_max = nslack_max;
    ws->N_max = N_max;

    ws->m = 0;
    ws->nslack = 0;
    ws->N = 0;

    ws->initialized = 1;

    /* double work arrays */
    ws->effective_rhs_data = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->x_B = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->lambda = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->d_dir = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->cB = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->w = (double *) c_ptr;
    c_ptr += ws->nv * sizeof(double);
    ws->lower = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->upper = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->nonbasic_value = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->singleton_value = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->feas_tol = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->dual_cost = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->dual_d = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->col_norm_sq = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->cert_q_scale = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->row_norm = (double *) c_ptr;
    c_ptr += m_max * sizeof(double);
    ws->lo_eff = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->devex_w = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->up_eff = (double *) c_ptr;
    c_ptr += N_max * sizeof(double);
    ws->has_lower = (unsigned char *) c_ptr;
    c_ptr += N_max * sizeof(unsigned char);
    ws->has_upper = (unsigned char *) c_ptr;
    c_ptr += N_max * sizeof(unsigned char);
    ws->status = (unsigned char *) c_ptr;
    c_ptr += N_max * sizeof(unsigned char);

    /* Align before the integer bookkeeping arrays. */
    align_char_to(8, &c_ptr);

    ws->basis_idx = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->lambda_nz = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->sing_nz = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->lambda_nnz = -1;
    ws->nonbasis_idx = (int *) c_ptr;
    c_ptr += N_max * sizeof(int);
    ws->ipiv = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->row_perm = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->singleton_row = (int *) c_ptr;
    c_ptr += N_max * sizeof(int);
    ws->singleton_basis_pos = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->singleton_rows = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->reduced_basis_pos = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->reduced_rows = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->row_to_singleton = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->basis_pos_to_reduced = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);
    ws->reduced_var = (int *) c_ptr;
    c_ptr += ws->nv * sizeof(int);
    ws->var_to_reduced = (int *) c_ptr;
    c_ptr += N_max * sizeof(int);
    ws->row_to_reduced = (int *) c_ptr;
    c_ptr += m_max * sizeof(int);

    /* align to 64 bytes before BLASFEO raw memory blocks */
    align_char_to(64, &c_ptr);

    /* BLASFEO matrices/vectors */
    blasfeo_create_dmat(m_max, dims->nv, &ws->Aphi, c_ptr);
    c_ptr += blasfeo_memsize_dmat(m_max, dims->nv);

    blasfeo_create_dmat(m_max, m_max, &ws->B_mat, c_ptr);
    c_ptr += blasfeo_memsize_dmat(m_max, m_max);

    blasfeo_create_dmat(m_max, m_max, &ws->LU_mat, c_ptr);
    c_ptr += blasfeo_memsize_dmat(m_max, m_max);

    blasfeo_create_dmat(m_max, dims->nv, &ws->G_mat, c_ptr);
    c_ptr += blasfeo_memsize_dmat(m_max, dims->nv);

    blasfeo_create_dmat(dims->nv, dims->nv, &ws->L_red, c_ptr);
    c_ptr += blasfeo_memsize_dmat(dims->nv, dims->nv);

    blasfeo_create_dmat(dims->nv, dims->nv, &ws->U_red, c_ptr);
    c_ptr += blasfeo_memsize_dmat(dims->nv, dims->nv);

    blasfeo_create_dvec(m_max, &ws->rhs, c_ptr);
    c_ptr += blasfeo_memsize_dvec(m_max);

    blasfeo_create_dvec(N_max, &ws->c, c_ptr);
    c_ptr += blasfeo_memsize_dvec(N_max);

    blasfeo_create_dvec(N_max, &ws->rc, c_ptr);
    c_ptr += blasfeo_memsize_dvec(N_max);

    blasfeo_create_dvec(m_max, &ws->tmp_m_1, c_ptr);
    c_ptr += blasfeo_memsize_dvec(m_max);

    blasfeo_create_dvec(m_max, &ws->tmp_m_2, c_ptr);
    c_ptr += blasfeo_memsize_dvec(m_max);

    blasfeo_create_dvec(m_max, &ws->lu_work, c_ptr);
    c_ptr += blasfeo_memsize_dvec(m_max);

    align_char_to(8, &c_ptr);

    *ws_out = ws;
    return c_ptr;
}

acados_size_t dense_lp_feasibility_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    (void) config_;
    (void) opts_;

    dense_qp_dims *dims = dims_;
    acados_size_t size = 0;

    size += sizeof(dense_lp_feasibility_memory);

    // no extra struct d_simplex_workspace here if already included in helper;
    // choose one convention and keep it consistent.
    size += dense_lp_feasibility_workspace_calculate_size_internal(dims);

    size += 8;
    make_int_multiple_of(8, &size);
    return size;
}

void *dense_lp_feasibility_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    (void) config_;
    (void) opts_;

    dense_qp_dims *dims = dims_;

    char *c_ptr = (char *) raw_memory;
    align_char_to(8, &c_ptr);

    dense_lp_feasibility_memory *mem = (dense_lp_feasibility_memory *) c_ptr;
    c_ptr += sizeof(dense_lp_feasibility_memory);

    align_char_to(8, &c_ptr);

    /* assign the private simplex workspace into the remaining memory */
    struct d_simplex_workspace *ws = NULL;
    c_ptr = dense_lp_feasibility_workspace_assign_internal(dims, c_ptr, &ws);
    mem->simplex_workspace = ws;

    mem->time_lp_solver_call = 0.0;
    mem->time_lp_setup = 0.0;
    mem->iter = 0;
    mem->status = 0;

    assert((char *) raw_memory + dense_lp_feasibility_memory_calculate_size(config_, dims_, opts_) >= c_ptr);

    return (void *) mem;
}

void dense_lp_feasibility_memory_get(void *config_, void *mem_, const char *field, void *value)
{
    (void) config_;

    dense_lp_feasibility_memory *mem = mem_;

    if (!strcmp(field, "time_lp_solver_call"))
    {
        *(double *) value = mem->time_lp_solver_call;
    }
    else if (!strcmp(field, "time_lp_setup"))
    {
        *(double *) value = mem->time_lp_setup;
    }
    else if (!strcmp(field, "iter"))
    {
        *(int *) value = mem->iter;
    }
    else if (!strcmp(field, "status"))
    {
        *(int *) value = mem->status;
    }
    else
    {
        printf("\nerror: dense_lp_feasibility_memory_get: field %s not available\n", field);
        exit(1);
    }
}

/************************************************
* misc
************************************************/

acados_size_t dense_lp_feasibility_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    (void) config_;
    (void) dims_;
    (void) opts_;

    /* no separate work object needed for now */
    return 0;
}

void dense_lp_feasibility_memory_reset(void *config, void *qp_in, void *qp_out, void *opts, void *mem, void *work)
{
    (void) config;
    (void) qp_in;
    (void) qp_out;
    (void) opts;
    (void) work;

    dense_lp_feasibility_memory *memory = mem;

    memory->time_lp_solver_call = 0.0;
    memory->time_lp_setup = 0.0;
    memory->iter = 0;
    memory->status = 0;
}

void dense_lp_feasibility_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_,
                                     const char *field, int stage, void *value, int size1, int size2)
{
    (void) config_;
    (void) qp_in_;
    (void) qp_out_;
    (void) opts_;
    (void) stage;
    (void) size1;
    (void) size2;

    dense_lp_feasibility_memory *mem = mem_;

    if (!strcmp(field, "status"))
    {
        *(int *) value = mem->status;
    }
    else if (!strcmp(field, "iter"))
    {
        *(int *) value = mem->iter;
    }
    else if (!strcmp(field, "time_lp_solver_call"))
    {
        *(double *) value = mem->time_lp_solver_call;
    }
    else if (!strcmp(field, "time_lp_setup"))
    {
        *(double *) value = mem->time_lp_setup;
    }
    else
    {
        printf("\nerror: dense_lp_feasibility_solver_get: field %s not available\n", field);
        exit(1);
    }
}

void dense_lp_feasibility_terminate(void *config_, void *mem_, void *work_)
{
    (void) config_;
    (void) mem_;
    (void) work_;
    return;
}

void dense_lp_feasibility_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->dims_set = &dense_qp_dims_set;

    config->opts_calculate_size = &dense_lp_feasibility_opts_calculate_size;
    config->opts_assign = &dense_lp_feasibility_opts_assign;
    config->opts_initialize_default = &dense_lp_feasibility_opts_initialize_default;
    config->opts_update = &dense_lp_feasibility_opts_update;
    config->opts_set = &dense_lp_feasibility_opts_set;
    config->opts_get = &dense_lp_feasibility_opts_get;

    config->memory_calculate_size = &dense_lp_feasibility_memory_calculate_size;
    config->memory_assign = &dense_lp_feasibility_memory_assign;
    config->memory_get = &dense_lp_feasibility_memory_get;

    config->workspace_calculate_size = &dense_lp_feasibility_workspace_calculate_size;
    config->evaluate = &dense_lp_feasibility;

    /* no sensitivities for this solver */
    config->eval_forw_sens = NULL;
    config->eval_adj_sens = NULL;

    config->memory_reset = &dense_lp_feasibility_memory_reset;
    config->solver_get = &dense_lp_feasibility_solver_get;
    config->terminate = &dense_lp_feasibility_terminate;
}


/* Copy the row in ws->rc (length nv) into Aphi row `row`, scaled to unit
 * norm if SIMPLEX_ROW_NORMALIZATION, accumulate the structural column norms
 * and set up the auxiliary variable of the row. Returns the row scale. */
static double d_simplex_setup_insert_row(struct d_simplex_workspace *ws, int nv, int row)
{
    double *x = ws->rc.pa;
    double norm_squared = 0.0;
    for (int j = 0; j < nv; ++j)
        norm_squared += x[j] * x[j];
    double row_scale = SIMPLEX_ROW_NORMALIZATION && norm_squared > 0.0 ?
        1.0 / sqrt(norm_squared) : 1.0;
    blasfeo_drowin(nv, row_scale, &ws->rc, 0, &ws->Aphi, row, 0);
    double s2 = row_scale * row_scale;
    for (int j = 0; j < nv; ++j)
        ws->col_norm_sq[j] += s2 * x[j] * x[j];

    int aux = nv + row;
    ws->row_norm[row] = row_scale * sqrt(norm_squared);
    ws->feas_tol[aux] = ws->row_norm[row] > 0.0 ?
        ws->tol_infeas * ws->row_norm[row] : ws->tol_infeas;
    ws->singleton_row[aux] = row;
    ws->singleton_value[aux] = -1.0;
    ws->nonbasic_value[aux] = 0.0;
    ws->status[aux] = D_SIMPLEX_BASIC;
    ws->basis_idx[row] = aux;
    return row_scale;
}

/*
 * Build the Phase-I LP from the dense QP. The QP data is assumed to be
 * consistent (checked by acados): only the bound masks are read, soft
 * constraints (idxs_rev >= 0) are dropped. Rows: equalities A x - s = b with
 * s fixed at 0, then one row c_j^T x - s_j = 0 with lg_j <= s_j <= ug_j per
 * hard generic constraint with at least one active side.
 */
static int d_simplex_workspace_update_from_qp(struct d_dense_qp *qp, struct d_simplex_workspace *ws)
{
    const int nv = qp->dim->nv;
    const int ne = qp->dim->ne;
    const int ng = qp->dim->ng;
    const int nb = qp->dim->nb;
    const double *d = qp->d->pa;
    const double *d_mask = qp->d_mask->pa;
    const int *idxs_rev = qp->idxs_rev;
    /* hpipm layout of d: [lb, lg, -ub, -ug] */
    const double *lbw = d, *lbg = d + nb, *ubw = d + nb + ng, *ubg = d + 2 * nb + ng;
    const double *lbw_m = d_mask, *lbg_m = d_mask + nb;
    const double *ubw_m = d_mask + nb + ng, *ubg_m = d_mask + 2 * nb + ng;

    /* structural variables */
    for (int j = 0; j < nv; ++j)
    {
        ws->has_lower[j] = 0;
        ws->has_upper[j] = 0;
        ws->lower[j] = 0.0;
        ws->upper[j] = 0.0;
        ws->singleton_row[j] = -1;
        ws->singleton_value[j] = 0.0;
        ws->feas_tol[j] = ws->tol_infeas;
        ws->col_norm_sq[j] = 0.0;
    }
    for (int k = 0; k < nb; ++k)
    {
        if (idxs_rev[k] >= 0)
            continue;
        const int idx = qp->idxb[k];
        if (lbw_m[k] != 0.0)
        {
            ws->has_lower[idx] = 1;
            ws->lower[idx] = lbw[k];
        }
        if (ubw_m[k] != 0.0)
        {
            ws->has_upper[idx] = 1;
            ws->upper[idx] = -ubw[k];
        }
        if (ws->has_lower[idx] && ws->has_upper[idx] && ws->lower[idx] > ws->upper[idx])
            return ACADOS_INFEASIBLE;
    }
    for (int j = 0; j < nv; ++j)
    {
        ws->nonbasis_idx[j] = j;
        if (ws->has_lower[j] && ws->has_upper[j] && ws->lower[j] == ws->upper[j])
        {
            ws->nonbasic_value[j] = ws->lower[j];
            ws->status[j] = D_SIMPLEX_FIXED_NONBASIC;
        }
        else if (!ws->has_lower[j] && !ws->has_upper[j])
        {
            ws->nonbasic_value[j] = 0.0;
            ws->status[j] = D_SIMPLEX_FREE_NONBASIC;
        }
        else if (SIMPLEX_SUPERBASIC_START &&
                 (!ws->has_lower[j] || ws->lower[j] < 0.0) &&
                 (!ws->has_upper[j] || ws->upper[j] > 0.0))
        {
            ws->nonbasic_value[j] = 0.0;
            ws->status[j] = D_SIMPLEX_SUPERBASIC;
        }
        else if (ws->has_lower[j] && (!SIMPLEX_SUPERBASIC_START || ws->lower[j] >= 0.0))
        {
            ws->nonbasic_value[j] = ws->lower[j];
            ws->status[j] = D_SIMPLEX_AT_LOWER;
        }
        else
        {
            ws->nonbasic_value[j] = ws->upper[j];
            ws->status[j] = D_SIMPLEX_AT_UPPER;
        }
    }

    /* rows and auxiliary variables; ws->rc is scratch during setup */
    int row = 0;
    for (int eq = 0; eq < ne; ++eq, ++row)
    {
        blasfeo_drowex(nv, 1.0, qp->A, eq, 0, &ws->rc, 0);
        double row_scale = d_simplex_setup_insert_row(ws, nv, row);
        int aux = nv + row;
        BLASFEO_DVECEL(&ws->rhs, row) = row_scale * BLASFEO_DVECEL(qp->b, eq);
        ws->has_lower[aux] = 1;
        ws->has_upper[aux] = 1;
        ws->lower[aux] = 0.0;
        ws->upper[aux] = 0.0;
    }
    for (int j = 0; j < ng; ++j)
    {
        int has_lb = lbg_m[j] != 0.0;
        int has_ub = ubg_m[j] != 0.0;
        if (idxs_rev[nb + j] >= 0 || (!has_lb && !has_ub))
            continue;
        if (has_lb && has_ub && lbg[j] > -ubg[j])
            return ACADOS_INFEASIBLE;
        blasfeo_dcolex(nv, qp->Ct, 0, j, &ws->rc, 0);
        double row_scale = d_simplex_setup_insert_row(ws, nv, row);
        int aux = nv + row;
        BLASFEO_DVECEL(&ws->rhs, row) = 0.0;
        ws->has_lower[aux] = has_lb;
        ws->has_upper[aux] = has_ub;
        ws->lower[aux] = has_lb ? row_scale * lbg[j] : 0.0;
        ws->upper[aux] = has_ub ? -row_scale * ubg[j] : 0.0;
        ++row;
    }

    const int N = nv + row;
    for (int j = nv; j < N; ++j)
        ws->nonbasis_idx[j] = -1;
    for (int j = 0; j < N; ++j)
    {
        ws->lo_eff[j] = ws->has_lower[j] ? ws->lower[j] : -INFINITY;
        ws->up_eff[j] = ws->has_upper[j] ? ws->upper[j] : INFINITY;
    }

    ws->nv = nv;
    ws->ne = ne;
    ws->ng = ng;
    ws->nb = nb;
    ws->m = row;
    ws->nslack = row;
    ws->nineq_slack = row - ne;
    ws->neq_residual = ne;
    ws->N = N;
    ws->iter = 0;
    ws->reduced_active = 0;
    ws->reduced_rep_valid = 0;
    ws->reduced_lu_valid = 0;
    blasfeo_dvecse(N, 0.0, &ws->c, 0);
    return ACADOS_SUCCESS;
}

/*
 * Acceptance test of the primal point ws->w: box bounds within tol and each
 * row value Aphi_i w - rhs_i within the bounds of its auxiliary variable,
 * up to tol * ||row_i||, tol = max(tol, tol_infeas). The comparisons are
 * written so that a NaN fails them.
 */
static int d_simplex_validate_original_solution(struct d_simplex_workspace *ws)
{
    const int nv = ws->nv, m = ws->m;
    const double tol = fmax(ws->tol, ws->tol_infeas);

    for (int j = 0; j < nv; ++j)
    {
        double x = ws->w[j];
        if (!(x == x) ||
            (ws->has_lower[j] && x < ws->lower[j] - tol) ||
            (ws->has_upper[j] && x > ws->upper[j] + tol))
            return 0;
    }
    if (m == 0)
        return 1;

    blasfeo_pack_dvec(nv, ws->w, 1, &ws->rc, 0);
    blasfeo_dgemv_n(m, nv, 1.0, &ws->Aphi, 0, 0, &ws->rc, 0,
                    -1.0, &ws->rhs, 0, &ws->tmp_m_1, 0);
    for (int i = 0; i < m; ++i)
    {
        int aux = nv + i;
        double value = BLASFEO_DVECEL(&ws->tmp_m_1, i);
        double row_tol = tol * (ws->row_norm[i] > 0.0 ? ws->row_norm[i] : 1.0);
        if (!(value == value) ||
            (ws->has_lower[aux] && value < ws->lower[aux] - row_tol) ||
            (ws->has_upper[aux] && value > ws->upper[aux] + row_tol))
            return 0;
    }
    return 1;
}

static void d_simplex_recover_soft_slacks(
    const struct d_dense_qp *qp,
    struct d_dense_qp_sol *qp_sol,
    const struct d_simplex_workspace *ws)
{
    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    if (ns == 0)
        return;

    blasfeo_dvecse(2 * ns, 0.0, qp_sol->v, nv);
    int slack_bound_offset = 2 * nb + 2 * ng;
    for (int slack = 0; slack < ns; slack++)
    {
        if (fabs(blasfeo_dvecex1(qp->d_mask, slack_bound_offset + slack)) > 0.5)
        {
            double lower_bound = blasfeo_dvecex1(qp->d, slack_bound_offset + slack);
            if (lower_bound > BLASFEO_DVECEL(qp_sol->v, nv + slack))
                BLASFEO_DVECEL(qp_sol->v, nv + slack) = lower_bound;
        }
        if (fabs(blasfeo_dvecex1(qp->d_mask, slack_bound_offset + ns + slack)) > 0.5)
        {
            double lower_bound = blasfeo_dvecex1(qp->d, slack_bound_offset + ns + slack);
            if (lower_bound > BLASFEO_DVECEL(qp_sol->v, nv + ns + slack))
                BLASFEO_DVECEL(qp_sol->v, nv + ns + slack) = lower_bound;
        }
    }

    for (int k = 0; k < nb; k++)
    {
        int slack = qp->idxs_rev[k];
        if (slack < 0)
            continue;

        double x = ws->w[qp->idxb[k]];
        if (fabs(blasfeo_dvecex1(qp->d_mask, k)) > 0.5)
        {
            double violation = fmax(0.0, blasfeo_dvecex1(qp->d, k) - x);
            double *slack_lower = &BLASFEO_DVECEL(qp_sol->v, nv + slack);
            if (violation > *slack_lower)
                *slack_lower = violation;
        }
        if (fabs(blasfeo_dvecex1(qp->d_mask, nb + ng + k)) > 0.5)
        {
            double upper = -blasfeo_dvecex1(qp->d, nb + ng + k);
            double violation = fmax(0.0, x - upper);
            double *slack_upper = &BLASFEO_DVECEL(qp_sol->v, nv + ns + slack);
            if (violation > *slack_upper)
                *slack_upper = violation;
        }
    }

    for (int j = 0; j < ng; j++)
    {
        int slack = qp->idxs_rev[nb + j];
        if (slack < 0)
            continue;

        double value = 0.0;
        for (int i = 0; i < nv; i++)
            value += blasfeo_dgeex1(qp->Ct, i, j) * ws->w[i];

        if (fabs(blasfeo_dvecex1(qp->d_mask, nb + j)) > 0.5)
        {
            double violation = fmax(0.0, blasfeo_dvecex1(qp->d, nb + j) - value);
            double *slack_lower = &BLASFEO_DVECEL(qp_sol->v, nv + slack);
            if (violation > *slack_lower)
                *slack_lower = violation;
        }
        if (fabs(blasfeo_dvecex1(qp->d_mask, 2 * nb + ng + j)) > 0.5)
        {
            double upper = -blasfeo_dvecex1(qp->d, 2 * nb + ng + j);
            double violation = fmax(0.0, value - upper);
            double *slack_upper = &BLASFEO_DVECEL(qp_sol->v, nv + ns + slack);
            if (violation > *slack_upper)
                *slack_upper = violation;
        }
    }
}

/*
 * Map every basic variable either to its position in the reduced basis or to
 * its singleton row. Returns an error on an inconsistent map or a NaN.
 */
static int d_simplex_refresh_reduced_basis_positions(
    struct d_simplex_workspace *ws,
    double tol_pivot)
{
    int m = ws->m;
    for (int pos = 0; pos < m; ++pos)
        ws->basis_pos_to_reduced[pos] = -1;
    for (int j = 0; j < ws->nreduced; ++j)
        ws->reduced_basis_pos[j] = -1;
    for (int p = 0; p < ws->nsingleton; ++p)
        ws->singleton_basis_pos[p] = -1;

    for (int pos = 0; pos < m; ++pos)
    {
        int var = ws->basis_idx[pos];
        int singleton_row = ws->singleton_row[var];
        if (singleton_row >= 0)
        {
            if (ws->row_to_singleton[singleton_row] < 0)
                return ACADOS_QP_FAILURE;
            int singleton_pos = ws->row_to_singleton[singleton_row];
            if (singleton_pos >= ws->nsingleton || ws->singleton_basis_pos[singleton_pos] >= 0)
                return ACADOS_QP_FAILURE;
            ws->singleton_basis_pos[singleton_pos] = pos;
        }
        else
        {
            int reduced_pos = ws->var_to_reduced[var];
            if (reduced_pos < 0 || reduced_pos >= ws->nreduced ||
                ws->reduced_basis_pos[reduced_pos] >= 0)
                return ACADOS_QP_FAILURE;
            ws->reduced_basis_pos[reduced_pos] = pos;
            ws->basis_pos_to_reduced[pos] = reduced_pos;
        }
    }

    for (int p = 0; p < ws->nsingleton; ++p)
    {
        if (ws->singleton_basis_pos[p] < 0)
            return ACADOS_QP_FAILURE;
    }
    for (int j = 0; j < ws->nreduced; ++j)
    {
        if (ws->reduced_basis_pos[j] < 0)
            return ACADOS_QP_FAILURE;
    }

    return ACADOS_SUCCESS;
}

/* Entry (row, var) of the Phase-I constraint matrix [A_struct, aux]. Only
 * the structural block is stored in Aphi; auxiliary column var >= nv is the
 * singleton singleton_value[var] * e_{singleton_row[var]}. */
static inline double d_simplex_A(const struct d_simplex_workspace *ws, int row, int var)
{
    if (var < ws->nv)
        return BLASFEO_DMATEL(&ws->Aphi, row, var);
    return row == ws->singleton_row[var] ? ws->singleton_value[var] : 0.0;
}

/* Column var of the Phase-I constraint matrix into the first m entries of x. */
static void d_simplex_column_extract(struct d_simplex_workspace *ws, int var,
                                     struct blasfeo_dvec *x)
{
    if (var < ws->nv)
    {
        blasfeo_dcolex(ws->m, &ws->Aphi, 0, var, x, 0);
        return;
    }
    blasfeo_dvecse(ws->m, 0.0, x, 0);
    BLASFEO_DVECEL(x, ws->singleton_row[var]) = ws->singleton_value[var];
}

/* Entry (i, j) of the transposed reduced basis T = H^T. */
static inline double d_simplex_reduced_T(const struct d_simplex_workspace *ws, int i, int j)
{
    return BLASFEO_DMATEL(&ws->Aphi, ws->reduced_rows[j], ws->reduced_var[i]);
}

#if SIMPLEX_DEBUG_LU_UPDATE
/* max |(P*T - L*U)_ij| / max(1, max |T_ij|) */
static double d_simplex_reduced_lu_residual(const struct d_simplex_workspace *ws)
{
    int k = ws->nreduced;
    double err = 0.0, scale = 1.0;
    for (int i = 0; i < k; ++i)
    {
        for (int j = 0; j < k; ++j)
        {
            double lu = 0.0;
            for (int l = 0; l <= (i < j ? i : j); ++l)
                lu += BLASFEO_DMATEL(&ws->L_red, i, l) * BLASFEO_DMATEL(&ws->U_red, l, j);
            double t = d_simplex_reduced_T(ws, ws->row_perm[i], j);
            err = fmax(err, fabs(lu - t));
            scale = fmax(scale, fabs(t));
        }
    }
    return err / scale;
}
#endif

/*
 * Fresh LU factorization with partial pivoting of the transposed reduced
 * basis T = H^T (gathered from Aphi through the reduced maps), split into
 * the separate factors L_red (unit lower) and U_red (upper) that the
 * Fletcher-Matthews update operates on.
 */
static int d_simplex_refactor_reduced_H(struct d_simplex_workspace *ws, double tol_pivot)
{
    int k = ws->nreduced;
#if SIMPLEX_PROFILING
    acados_timer timer;
#endif

    ws->reduced_lu_valid = 0;
    ws->lu_updates_since_refactor = 0;
    ws->reduced_refactor_count++;
    if (k > 0)
    {
#if SIMPLEX_PROFILING
        acados_tic(&timer);
#endif
        for (int j = 0; j < k; ++j)
        {
            int row = ws->reduced_rows[j];
            for (int i = 0; i < k; ++i)
            {
                BLASFEO_DMATEL(&ws->B_mat, i, j) =
                    BLASFEO_DMATEL(&ws->Aphi, row, ws->reduced_var[i]);
            }
        }
        ws->reduced_H_scalar_reads += k * k;
#if SIMPLEX_PROFILING
        ws->time_reduced_H_build += acados_toc(&timer);
        acados_tic(&timer);
#endif
        blasfeo_dgetrf_rp(k, k, &ws->B_mat, 0, 0, &ws->LU_mat, 0, 0, ws->ipiv);
#if SIMPLEX_PROFILING
        ws->time_reduced_factorization_only += acados_toc(&timer);
#endif

        for (int i = 0; i < k; ++i)
            ws->row_perm[i] = i;
        for (int i = 0; i < k; ++i)
        {
            int pivot = ws->ipiv[i];
            /* negated test: also rejects a NaN pivot */
            if (!(fabs(BLASFEO_DMATEL(&ws->LU_mat, i, i)) > tol_pivot))
                return ACADOS_NAN_DETECTED;
            int tmp = ws->row_perm[i];
            ws->row_perm[i] = ws->row_perm[pivot];
            ws->row_perm[pivot] = tmp;
        }

        for (int j = 0; j < k; ++j)
        {
            for (int i = 0; i < k; ++i)
            {
                double value = BLASFEO_DMATEL(&ws->LU_mat, i, j);
                BLASFEO_DMATEL(&ws->L_red, i, j) = i > j ? value : (i == j ? 1.0 : 0.0);
                BLASFEO_DMATEL(&ws->U_red, i, j) = i <= j ? value : 0.0;
            }
        }
        ws->L_red.use_dA = 0;
        ws->U_red.use_dA = 0;
    }

    ws->reduced_lu_valid = 1;
    return ACADOS_SUCCESS;
}

/*
 * Delete column q of the nrows x ncols block of U: shift columns
 * q+1..ncols-1 one position left and clear column ncols-1. In both BLASFEO
 * layouts a block of consecutive columns restricted to `col_stride` rows (one
 * panel in panel-major, all rows in column-major) is contiguous, so the shift
 * is one memmove per row block.
 */
static void d_simplex_lu_delete_column(int nrows, int ncols, int q, struct blasfeo_dmat *U)
{
    if (q < ncols - 1)
    {
        const ptrdiff_t col_stride = &BLASFEO_DMATEL(U, 0, 1) - &BLASFEO_DMATEL(U, 0, 0);
        const size_t block = (size_t) (ncols - 1 - q) * (size_t) col_stride * sizeof(double);
        for (int r0 = 0; r0 < nrows; r0 += (int) col_stride)
            memmove(&BLASFEO_DMATEL(U, r0, q), &BLASFEO_DMATEL(U, r0, q + 1), block);
    }
    blasfeo_dgese(nrows, 1, 0.0, U, 0, ncols - 1);
}

/*
 * Fletcher-Matthews restoration of perm * A = L * U after a column deletion
 * (R. Fletcher, F.J. Matthews, "Stable modification of explicit LU factors
 * for simplex updates", Math. Program. 30, 1984; two-case elimination and
 * notation u1ii, u2ii, l21 follow the paper):
 * U is upper Hessenberg in columns q..m-2 (entries U(ii+1, ii)), L is unit
 * lower triangular. Each subdiagonal entry is eliminated by a row operation,
 * with an interchange of rows ii, ii+1 (of perm, U and the leading part of L)
 * when that bounds the multiplier by 1. U has m rows and at least m-1
 * columns (columns beyond m-1 are carried along). O(m^2).
 */
static int d_simplex_fm_restore(
    int m, int ncols, int q,
    struct blasfeo_dmat *L,
    struct blasfeo_dmat *U,
    int *perm,
    struct blasfeo_dvec *tmp1,
    struct blasfeo_dvec *tmp2,
    double tol_pivot)
{
    for (int ii = q; ii < m - 1; ii++)
    {
        double u1ii = BLASFEO_DMATEL(U, ii, ii);
        double u2ii = BLASFEO_DMATEL(U, ii + 1, ii);
        double l21 = BLASFEO_DMATEL(L, ii + 1, ii);
        double delta = u2ii + l21 * u1ii;
        int tail = ncols - ii - 1;

        if (fabs(u1ii) > fabs(delta))
        {
            /* Fletcher-Matthews case 1: eliminate without interchange. */
            if (!(fabs(u1ii) > tol_pivot))
                return ACADOS_NAN_DETECTED;
            double r = u2ii / u1ii;

            if (tail > 0)
            {
                blasfeo_drowex(tail, 1.0, U, ii, ii + 1, tmp1, 0);
                blasfeo_drowad(tail, -r, tmp1, 0, U, ii + 1, ii + 1);
            }
            int below = m - ii - 1;
            if (below > 0)
            {
                blasfeo_dcolex(below, L, ii + 1, ii + 1, tmp1, 0);
                blasfeo_dcolad(below, r, tmp1, 0, L, ii + 1, ii);
            }
            BLASFEO_DMATEL(U, ii + 1, ii) = 0.0;
            BLASFEO_DMATEL(L, ii, ii) = 1.0;
            BLASFEO_DMATEL(L, ii + 1, ii + 1) = 1.0;
        }
        else
        {
            /* Fletcher-Matthews case 2: interchange rows ii and ii+1. */
            if (!(fabs(delta) > tol_pivot))
                return ACADOS_NAN_DETECTED;
            double alpha = u1ii / delta;

            if (ii > 0)
                blasfeo_drowsw(ii, L, ii, 0, L, ii + 1, 0);

            int tmp_perm = perm[ii];
            perm[ii] = perm[ii + 1];
            perm[ii + 1] = tmp_perm;

            /* new_u1 = old_u2 + l21*old_u1, new_u2 = old_u1 - alpha*new_u1 */
            if (tail > 0)
            {
                blasfeo_drowex(tail, 1.0, U, ii, ii + 1, tmp1, 0);
                blasfeo_drowex(tail, 1.0, U, ii + 1, ii + 1, tmp2, 0);
                blasfeo_daxpy(tail, l21, tmp1, 0, tmp2, 0, tmp2, 0);
                blasfeo_daxpy(tail, -alpha, tmp2, 0, tmp1, 0, tmp1, 0);
                blasfeo_drowin(tail, 1.0, tmp2, 0, U, ii, ii + 1);
                blasfeo_drowin(tail, 1.0, tmp1, 0, U, ii + 1, ii + 1);
            }
            BLASFEO_DMATEL(U, ii, ii) = delta;
            BLASFEO_DMATEL(U, ii + 1, ii) = 0.0;

            /* new_l2 = old_l1 - l21*old_l2, new_l1 = old_l2 + alpha*new_l2 */
            int below = m - ii - 2;
            if (below > 0)
            {
                blasfeo_dcolex(below, L, ii + 2, ii, tmp1, 0);
                blasfeo_dcolex(below, L, ii + 2, ii + 1, tmp2, 0);
                blasfeo_daxpy(below, -l21, tmp2, 0, tmp1, 0, tmp1, 0);
                blasfeo_daxpy(below, alpha, tmp1, 0, tmp2, 0, tmp2, 0);
                blasfeo_dcolin(below, tmp2, 0, L, ii + 2, ii);
                blasfeo_dcolin(below, tmp1, 0, L, ii + 2, ii + 1);
            }
            BLASFEO_DMATEL(L, ii, ii) = 1.0;
            BLASFEO_DMATEL(L, ii + 1, ii) = alpha;
            BLASFEO_DMATEL(L, ii + 1, ii + 1) = 1.0;
        }
    }
    L->use_dA = 0;
    U->use_dA = 0;
    return ACADOS_SUCCESS;
}

/* Negated test, also rejects a NaN pivot. */
static int d_simplex_lu_diag_ok(int m, const struct blasfeo_dmat *U, double tol_pivot)
{
    for (int i = 0; i < m; i++)
        if (!(fabs(BLASFEO_DMATEL(U, i, i)) > tol_pivot))
            return 0;
    return 1;
}

/*
 * Fletcher-Matthews update of perm * A = L * U for a column replacement in
 * delete-and-append form:
 *
 *     A_new = [A(:,0:q-1), A(:,q+1:m-1), a_enter].
 *
 * U(:,q) is deleted, the following columns are shifted left, the upper
 * Hessenberg matrix is restored to triangular form (d_simplex_fm_restore),
 * and the new column is L^{-1} perm a_enter. O(m^2). work, tmp1 and tmp2
 * must hold at least m entries. On failure, L, U and perm are in an
 * undefined state and a refactorization is required.
 */
static int d_simplex_lu_delete_append_update(
    int m, int q,
    struct blasfeo_dmat *L,
    struct blasfeo_dmat *U,
    int *perm,
    struct blasfeo_dvec *a_enter,
    struct blasfeo_dvec *work,
    struct blasfeo_dvec *tmp1,
    struct blasfeo_dvec *tmp2,
    double tol_pivot)
{
    if (m <= 0 || q < 0 || q >= m)
        return ACADOS_NAN_DETECTED;

    d_simplex_lu_delete_column(m, m, q, U);
    if (d_simplex_fm_restore(m, m, q, L, U, perm, tmp1, tmp2, tol_pivot) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;

    /* last column of U: L^{-1} * perm * a_enter */
    blasfeo_dvecex_sp(m, 1.0, perm, a_enter, 0, work, 0);
    blasfeo_dtrsv_lnu(m, L, 0, 0, work, 0, work, 0);
    blasfeo_dcolin(m, work, 0, U, 0, m - 1);
    U->use_dA = 0;

    return d_simplex_lu_diag_ok(m, U, tol_pivot) ? ACADOS_SUCCESS : ACADOS_NAN_DETECTED;
}

/* TODO: implement a mirrored Fletcher-Matthews loop that works directly on
 * the transposed indices of L and U, to avoid the O(k^2) copies into and out
 * of the transposed view, which dominate the update cost for large k. */
/*
 * Transposed view of P T = L U (T is n x nc after a column deletion, U is
 * n x n, L is nc x nc of which the first n columns are used):
 *
 *     T^T P^T = U^T L^T = (U^T D^{-1}) (D L^T) = Lh Uh,  D = diag(U),
 *
 * with Lh (n x n) unit lower and Uh (n x nc) upper, the form required by the
 * Fletcher-Matthews routines: a row change of T is a column change of T^T.
 * Uses BLASFEO triangular transposes and column scalings (D L^T = (L D)^T),
 * so L is overwritten (scaled); only the triangles of Lh / Uh are written,
 * which is all that the Fletcher-Matthews routines and the triangular solves
 * read. tmp holds at least n entries.
 */
static void d_simplex_lu_to_transposed(
    int n, int nc,
    struct blasfeo_dmat *L, const struct blasfeo_dmat *U,
    struct blasfeo_dmat *Lh, struct blasfeo_dmat *Uh, struct blasfeo_dvec *tmp)
{
    if (n <= 0)
        return;
    blasfeo_dtrtr_u(n, (struct blasfeo_dmat *) U, 0, 0, Lh, 0, 0);   /* Lh = U^T */
    for (int j = 0; j < n; ++j)
    {
        double d = BLASFEO_DMATEL(U, j, j);
        blasfeo_dcolsc(n - j, 1.0 / d, Lh, j, j);                     /* U^T D^{-1} */
        blasfeo_dcolsc(nc - j, d, L, j, j);                           /* L D */
    }
    blasfeo_dtrtr_l(n, L, 0, 0, Uh, 0, 0);                            /* Uh = (L D)^T */
    if (nc > n)
    {
        /* rectangular case nc = n + 1: last row of L D is the last column */
        blasfeo_drowex(n, 1.0, L, n, 0, tmp, 0);
        blasfeo_dcolin(n, tmp, 0, Uh, 0, n);
    }
    Lh->use_dA = 0;
    Uh->use_dA = 0;
    L->use_dA = 0;
}

/* Inverse of d_simplex_lu_to_transposed for square n x n factors:
 * perm_h T^T Q = Lh Uh  gives  Q^T (T perm_h^T) = L U with
 * L = Uh^T Dh^{-1}, U = Dh Lh^T = (Lh Dh)^T, Dh = diag(Uh).
 * Lh is overwritten (scaled). */
static void d_simplex_lu_from_transposed(
    int n,
    struct blasfeo_dmat *Lh, const struct blasfeo_dmat *Uh,
    struct blasfeo_dmat *L, struct blasfeo_dmat *U)
{
    if (n <= 0)
        return;
    blasfeo_dtrtr_u(n, (struct blasfeo_dmat *) Uh, 0, 0, L, 0, 0);   /* L = Uh^T */
    for (int j = 0; j < n; ++j)
    {
        double d = BLASFEO_DMATEL(Uh, j, j);
        blasfeo_dcolsc(n - j, 1.0 / d, L, j, j);                      /* Uh^T Dh^{-1} */
        blasfeo_dcolsc(n - j, d, Lh, j, j);                           /* Lh Dh */
    }
    blasfeo_dtrtr_l(n, Lh, 0, 0, U, 0, 0);                            /* U = (Lh Dh)^T */
    L->use_dA = 0;
    U->use_dA = 0;
    Lh->use_dA = 0;
}

/*
 * Bordered extension of P*T = L*U after a structural variable entered the
 * reduced basis at position old_k together with a new reduced row at
 * position old_k (the reduced maps must already be updated):
 *
 *     T_new = [T  b; c^T d],  L_new = [L 0; l^T 1],  U_new = [U u; 0 s],
 *     L u = P b,  U^T l = c,  s = d - l^T u,  perm[old_k] = old_k.
 * O(k^2), no pivoting: a small Schur complement s triggers a refactorization.
 */
static int d_simplex_extend_reduced_lu(
    struct d_simplex_workspace *ws,
    int old_k,
    double tol_pivot)
{
    if (old_k < 0 || old_k >= ws->nreduced)
        return ACADOS_QP_FAILURE;

    struct blasfeo_dmat *L = &ws->L_red;
    struct blasfeo_dmat *U = &ws->U_red;

    for (int i = 0; i < old_k; ++i)
        BLASFEO_DVECEL(&ws->tmp_m_1, i) = d_simplex_reduced_T(ws, i, old_k);
    blasfeo_dvecex_sp(old_k, 1.0, ws->row_perm, &ws->tmp_m_1, 0, &ws->tmp_m_2, 0);
    if (old_k > 0)
        blasfeo_dtrsv_lnu(old_k, L, 0, 0, &ws->tmp_m_2, 0, &ws->tmp_m_2, 0);

    for (int j = 0; j < old_k; ++j)
        BLASFEO_DVECEL(&ws->tmp_m_1, j) = d_simplex_reduced_T(ws, old_k, j);
    if (old_k > 0)
        blasfeo_dtrsv_utn(old_k, U, 0, 0, &ws->tmp_m_1, 0, &ws->lu_work, 0);

    double schur = d_simplex_reduced_T(ws, old_k, old_k);
    for (int i = 0; i < old_k; ++i)
    {
        double u_i = BLASFEO_DVECEL(&ws->tmp_m_2, i);
        double l_i = BLASFEO_DVECEL(&ws->lu_work, i);
        schur -= l_i * u_i;
        BLASFEO_DMATEL(U, i, old_k) = u_i;
        BLASFEO_DMATEL(U, old_k, i) = 0.0;
        BLASFEO_DMATEL(L, old_k, i) = l_i;
        BLASFEO_DMATEL(L, i, old_k) = 0.0;
    }
    if (!(fabs(schur) > tol_pivot))
        return ACADOS_NAN_DETECTED;

    BLASFEO_DMATEL(U, old_k, old_k) = schur;
    BLASFEO_DMATEL(L, old_k, old_k) = 1.0;
    ws->row_perm[old_k] = old_k;
    L->use_dA = 0;
    U->use_dA = 0;
    return ACADOS_SUCCESS;
}

/*
 * Split the basis into basic slacks with a singleton row and reduced
 * (structural) variables, build the reduced maps and G_mat (singleton rows x
 * reduced variables), and factorize the reduced block. Returns an error on a
 * NaN or an inconsistent basis.
 */
static int d_simplex_build_reduced_factorization(struct d_simplex_workspace *ws, double tol_pivot)
{
    int m = ws->m;
    int ns = 0;
    int k = 0;
#if SIMPLEX_PROFILING
    acados_timer timer;
#endif

#if SIMPLEX_PROFILING
    acados_tic(&timer);
#endif
    for (int row = 0; row < m; ++row)
    {
        ws->row_to_singleton[row] = -1;
        ws->row_to_reduced[row] = -1;
    }
    for (int var = 0; var < ws->N_max; ++var)
        ws->var_to_reduced[var] = -1;

    for (int pos = 0; pos < m; ++pos)
    {
        int var = ws->basis_idx[pos];
        int singleton_row = ws->singleton_row[var];
        if (singleton_row >= 0)
        {
            ws->singleton_basis_pos[ns] = pos;
            ws->singleton_rows[ns] = singleton_row;
            ws->row_to_singleton[singleton_row] = ns++;
        }
        else
        {
            ws->reduced_var[k] = var;
            ws->var_to_reduced[var] = k;
            ws->reduced_basis_pos[k] = pos;
            k++;
        }
    }

    int nr = 0;
    for (int row = 0; row < m; ++row)
    {
        if (ws->row_to_singleton[row] < 0)
        {
            ws->reduced_rows[nr] = row;
            ws->row_to_reduced[row] = nr++;
        }
    }
    if (nr != k || ns + k != m)
        return ACADOS_NAN_DETECTED;

    ws->nsingleton = ns;
    ws->nreduced = k;
#if SIMPLEX_PROFILING
    ws->time_reduced_maps += acados_toc(&timer);
#endif

    if (d_simplex_refresh_reduced_basis_positions(ws, tol_pivot) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;

#if SIMPLEX_PROFILING
    acados_tic(&timer);
#endif
    for (int p = 0; p < ns; ++p)
    {
        int physical_row = ws->singleton_rows[p];
        for (int j = 0; j < k; ++j)
        {
            int structural_var = ws->reduced_var[j];
            BLASFEO_DMATEL(&ws->G_mat, p, j) =
                BLASFEO_DMATEL(&ws->Aphi, physical_row, structural_var);
            ws->reduced_G_scalar_reads++;
        }
    }
#if SIMPLEX_PROFILING
    ws->time_reduced_G_build += acados_toc(&timer);
#endif

    ws->reduced_rep_valid = 1;
    return d_simplex_refactor_reduced_H(ws, tol_pivot);
}

/*
 * Row replacement in P T = L U (pivot type 3: structural q replaces the
 * structural at reduced position p; the maps are already updated, so row p
 * of T is the new row t_j = Aphi(reduced_rows[j], q)). In the transposed
 * view this is a Fletcher-Matthews delete-append of column c0 (the position
 * of row p in P); its row interchanges reorder the reduced rows (columns of
 * T). O(k^2).
 */
static int d_simplex_lu_replace_row(struct d_simplex_workspace *ws, int p, int q, double tol_pivot)
{
    const int k = ws->nreduced;
    int *perm_h = ws->ipiv;      /* scratch */
    int *new_rows = ws->sing_nz; /* scratch */
    int c0 = -1;
    for (int i = 0; i < k; ++i)
        if (ws->row_perm[i] == p)
            c0 = i;
    if (c0 < 0)
        return ACADOS_QP_FAILURE;

    d_simplex_lu_to_transposed(k, k, &ws->L_red, &ws->U_red, &ws->B_mat, &ws->LU_mat, &ws->tmp_m_1);
    for (int j = 0; j < k; ++j)
    {
        perm_h[j] = j;
        BLASFEO_DVECEL(&ws->rc, j) = BLASFEO_DMATEL(&ws->Aphi, ws->reduced_rows[j], q);
    }
    if (d_simplex_lu_delete_append_update(k, c0, &ws->B_mat, &ws->LU_mat, perm_h, &ws->rc,
            &ws->lu_work, &ws->tmp_m_1, &ws->tmp_m_2, tol_pivot) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;
    d_simplex_lu_from_transposed(k, &ws->B_mat, &ws->LU_mat, &ws->L_red, &ws->U_red);

    /* rows of T in the new order: row p moves to the end */
    for (int i = c0; i < k - 1; ++i)
        ws->row_perm[i] = ws->row_perm[i + 1];
    ws->row_perm[k - 1] = p;
    /* columns of T (reduced rows) reordered by the interchanges */
    for (int j = 0; j < k; ++j)
        new_rows[j] = ws->reduced_rows[perm_h[j]];
    for (int j = 0; j < k; ++j)
    {
        ws->reduced_rows[j] = new_rows[j];
        ws->row_to_reduced[new_rows[j]] = j;
    }
    return d_simplex_lu_diag_ok(k, &ws->U_red, tol_pivot) ? ACADOS_SUCCESS : ACADOS_NAN_DETECTED;
}

/*
 * Shrink P T = L U by row p and column q (pivot type 2: the auxiliary of the
 * reduced row at column q enters, the structural at reduced position p
 * leaves). Called before the maps are updated. Column q is deleted with a
 * Fletcher-Matthews restoration of the normal factors (the last row of U
 * becomes zero), row p with one of the transposed factors. On success the
 * new row order of T is in row_perm[0:k-1] (still with the old position
 * numbers) and the new order of the remaining reduced-row positions in
 * col_order[0:k-1] (positions of the old reduced_rows with q removed).
 * O(k^2).
 */
static int d_simplex_lu_delete_row_column(struct d_simplex_workspace *ws, int p, int q,
                                          int *col_order, double tol_pivot)
{
    const int k = ws->nreduced;
    int *perm_h = ws->ipiv; /* scratch */
    if (k <= 1)
        return ACADOS_SUCCESS;

    /* 1) delete column q of T: perm T(:, ~q) = L(:, 0:k-2) U(0:k-2, 0:k-2) */
    d_simplex_lu_delete_column(k, k, q, &ws->U_red);
    if (d_simplex_fm_restore(k, k - 1, q, &ws->L_red, &ws->U_red, ws->row_perm,
            &ws->tmp_m_1, &ws->tmp_m_2, tol_pivot) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;
    if (!d_simplex_lu_diag_ok(k - 1, &ws->U_red, tol_pivot))
        return ACADOS_NAN_DETECTED;

    /* 2) delete row p: column c0 of the transposed factors */
    int c0 = -1;
    for (int i = 0; i < k; ++i)
        if (ws->row_perm[i] == p)
            c0 = i;
    if (c0 < 0)
        return ACADOS_QP_FAILURE;
    d_simplex_lu_to_transposed(k - 1, k, &ws->L_red, &ws->U_red, &ws->B_mat, &ws->LU_mat, &ws->tmp_m_1);
    d_simplex_lu_delete_column(k - 1, k, c0, &ws->LU_mat);
    for (int j = 0; j < k - 1; ++j)
        perm_h[j] = j;
    if (d_simplex_fm_restore(k - 1, k - 1, c0, &ws->B_mat, &ws->LU_mat, perm_h,
            &ws->tmp_m_1, &ws->tmp_m_2, tol_pivot) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;
    if (!d_simplex_lu_diag_ok(k - 1, &ws->LU_mat, tol_pivot))
        return ACADOS_NAN_DETECTED;
    d_simplex_lu_from_transposed(k - 1, &ws->B_mat, &ws->LU_mat, &ws->L_red, &ws->U_red);

    for (int i = c0; i < k - 1; ++i)
        ws->row_perm[i] = ws->row_perm[i + 1];
    /* remaining columns: old positions without q, reordered by perm_h */
    for (int j = 0; j < k - 1; ++j)
    {
        int pos = perm_h[j];
        col_order[j] = pos < q ? pos : pos + 1;
    }
    return ACADOS_SUCCESS;
}

/*
 * Basis-change bookkeeping of the reduced representation.
 */
/*
 * Update the reduced maps and the factors of T = H^T after a basis change.
 *
 *  type 0 (singleton enters / singleton leaves): the reduced row of the
 *         entering auxiliary is replaced by the row of the leaving one, i.e.
 *         a column of T is replaced -> Fletcher-Matthews update.
 *  type 1 (structural enters / singleton leaves): k grows by one ->
 *         bordered extension.
 *  type 3 (structural enters / structural leaves): a row of T is replaced
 *         -> d_simplex_lu_replace_row.
 *  type 2 (singleton enters / structural leaves): k shrinks by one, T loses
 *         a row and a column -> d_simplex_lu_delete_row_column.
 *
 * Returns ACADOS_SUCCESS if the maps are consistent; ws->reduced_lu_valid
 * tells whether the factors were updated or need a refactorization.
 */
static int d_simplex_update_reduced_representation(
    struct d_simplex_workspace *ws,
    int pivot_type,
    int entering_idx,
    int leaving_idx,
    double tol_pivot)
{
    int k = ws->nreduced;
    int ns = ws->nsingleton;
    int lu_was_valid = ws->reduced_lu_valid;
    int lu_status = ACADOS_QP_FAILURE;
#if SIMPLEX_PROFILING
    acados_timer timer;
    acados_tic(&timer);
#endif
    ws->reduced_lu_valid = 0;

    if (pivot_type == 0)
    {
        int entering_row = ws->singleton_row[entering_idx];
        int leaving_row = ws->singleton_row[leaving_idx];
        if (entering_row < 0 || leaving_row < 0 ||
            ws->row_to_reduced[entering_row] < 0 ||
            ws->row_to_singleton[leaving_row] < 0)
            return ACADOS_QP_FAILURE;

        int reduced_row_pos = ws->row_to_reduced[entering_row];
        int singleton_pos = ws->row_to_singleton[leaving_row];
        for (int j = 0; j < k; ++j)
            BLASFEO_DMATEL(&ws->G_mat, singleton_pos, j) =
                BLASFEO_DMATEL(&ws->Aphi, entering_row, ws->reduced_var[j]);

        /* delete-and-append order of the reduced rows (= columns of T) */
        for (int i = reduced_row_pos; i < k - 1; ++i)
        {
            ws->reduced_rows[i] = ws->reduced_rows[i + 1];
            ws->row_to_reduced[ws->reduced_rows[i]] = i;
        }
        ws->reduced_rows[k - 1] = leaving_row;
        ws->row_to_reduced[leaving_row] = k - 1;
        ws->row_to_reduced[entering_row] = -1;
        ws->row_to_singleton[entering_row] = singleton_pos;
        ws->singleton_rows[singleton_pos] = entering_row;
        ws->row_to_singleton[leaving_row] = -1;

        if (lu_was_valid)
        {
            /* new column of T: row of the leaving auxiliary */
            for (int i = 0; i < k; ++i)
                BLASFEO_DVECEL(&ws->rc, i) =
                    BLASFEO_DMATEL(&ws->Aphi, leaving_row, ws->reduced_var[i]);
            lu_status = d_simplex_lu_delete_append_update(
                k, reduced_row_pos, &ws->L_red, &ws->U_red, ws->row_perm,
                &ws->rc, &ws->lu_work, &ws->tmp_m_1, &ws->tmp_m_2, tol_pivot);
        }
    }
    else if (pivot_type == 1)
    {
        int leaving_row = ws->singleton_row[leaving_idx];
        int singleton_pos = leaving_row >= 0 ? ws->row_to_singleton[leaving_row] : -1;
        if (entering_idx >= ws->nv || leaving_row < 0 || singleton_pos < 0 ||
            k >= ws->nv || ws->var_to_reduced[entering_idx] >= 0)
            return ACADOS_QP_FAILURE;

        int last_singleton = ns - 1;
        if (singleton_pos != last_singleton)
        {
            for (int j = 0; j < k; ++j)
            {
                double value = BLASFEO_DMATEL(&ws->G_mat, last_singleton, j);
                BLASFEO_DMATEL(&ws->G_mat, singleton_pos, j) = value;
            }
            int moved_row = ws->singleton_rows[last_singleton];
            ws->singleton_rows[singleton_pos] = moved_row;
            ws->row_to_singleton[moved_row] = singleton_pos;
        }
        ws->row_to_singleton[leaving_row] = -1;
        ws->row_to_reduced[leaving_row] = k;
        ws->reduced_rows[k] = leaving_row;
        ws->reduced_var[k] = entering_idx;
        ws->var_to_reduced[entering_idx] = k;
        ws->nreduced = k + 1;
        ws->nsingleton = ns - 1;

        for (int p = 0; p < ws->nsingleton; ++p)
        {
            int row = ws->singleton_rows[p];
            BLASFEO_DMATEL(&ws->G_mat, p, k) =
                BLASFEO_DMATEL(&ws->Aphi, row, entering_idx);
        }

        if (lu_was_valid)
            lu_status = d_simplex_extend_reduced_lu(ws, k, tol_pivot);
    }
    else if (pivot_type == 3)
    {
        int reduced_col_pos = ws->var_to_reduced[leaving_idx];
        if (entering_idx >= ws->nv || reduced_col_pos < 0 ||
            ws->var_to_reduced[entering_idx] >= 0)
            return ACADOS_QP_FAILURE;

        for (int p = 0; p < ns; ++p)
        {
            int row = ws->singleton_rows[p];
            BLASFEO_DMATEL(&ws->G_mat, p, reduced_col_pos) =
                BLASFEO_DMATEL(&ws->Aphi, row, entering_idx);
        }
        ws->var_to_reduced[leaving_idx] = -1;
        ws->var_to_reduced[entering_idx] = reduced_col_pos;
        ws->reduced_var[reduced_col_pos] = entering_idx;
        /* row replacement in T */
        if (lu_was_valid)
            lu_status = d_simplex_lu_replace_row(ws, reduced_col_pos, entering_idx, tol_pivot);
    }
    else if (pivot_type == 2)
    {
        /* the auxiliary of reduced row e enters, the structural at reduced
         * position p leaves: T loses row p and column row_to_reduced[e] */
        int e = ws->singleton_row[entering_idx];
        int p = ws->var_to_reduced[leaving_idx];
        if (e < 0 || p < 0 || ws->row_to_reduced[e] < 0)
            return ACADOS_QP_FAILURE;
        int q = ws->row_to_reduced[e];
        int *col_order = ws->lambda_nz; /* scratch (lambda pattern is recomputed) */
        ws->lambda_nnz = -1;

        if (lu_was_valid)
            lu_status = d_simplex_lu_delete_row_column(ws, p, q, col_order, tol_pivot);
        if (lu_status != ACADOS_SUCCESS)
        {
            /* natural order; the factors are rebuilt from the maps */
            for (int j = 0; j < k - 1; ++j)
                col_order[j] = j < q ? j : j + 1;
        }

        /* reduced rows: drop e, apply the new column order */
        int *new_rows = ws->sing_nz; /* scratch */
        for (int j = 0; j < k - 1; ++j)
            new_rows[j] = ws->reduced_rows[col_order[j]];
        for (int j = 0; j < k - 1; ++j)
        {
            ws->reduced_rows[j] = new_rows[j];
            ws->row_to_reduced[new_rows[j]] = j;
        }
        ws->row_to_reduced[e] = -1;

        /* reduced variables: move the last one into position p */
        int last = k - 1;
        ws->var_to_reduced[leaving_idx] = -1;
        if (p != last)
        {
            int moved = ws->reduced_var[last];
            ws->reduced_var[p] = moved;
            ws->var_to_reduced[moved] = p;
            for (int r = 0; r < ns; ++r)
                BLASFEO_DMATEL(&ws->G_mat, r, p) = BLASFEO_DMATEL(&ws->G_mat, r, last);
            if (lu_status == ACADOS_SUCCESS)
                for (int i = 0; i < k - 1; ++i)
                    if (ws->row_perm[i] == last)
                        ws->row_perm[i] = p;
        }
        ws->nreduced = k - 1;

        /* row e becomes a singleton row: new row of G */
        for (int j = 0; j < k - 1; ++j)
            BLASFEO_DMATEL(&ws->G_mat, ns, j) =
                BLASFEO_DMATEL(&ws->Aphi, e, ws->reduced_var[j]);
        ws->singleton_rows[ns] = e;
        ws->row_to_singleton[e] = ns;
        ws->nsingleton = ns + 1;
    }
    else
    {
        return ACADOS_QP_FAILURE;
    }

    if (lu_status == ACADOS_SUCCESS)
    {
        ws->lu_update_count++;
        ws->reduced_lu_valid =
            ++ws->lu_updates_since_refactor < SIMPLEX_LU_UPDATE_MAX;
#if SIMPLEX_DEBUG_LU_UPDATE
        double residual = d_simplex_reduced_lu_residual(ws);
        if (residual > 1e-10)
            printf("d_simplex: LU update (type %d) residual %.3e\n", pivot_type, residual);
#endif
    }
    else if (lu_was_valid && pivot_type != 3)
    {
        ws->lu_update_fail_count++;
    }

    ws->G_mat.use_dA = 0;
#if SIMPLEX_PROFILING
    ws->time_reduced_incremental_update += acados_toc(&timer);
#endif
    return ACADOS_SUCCESS;
}

static void d_simplex_reduced_primal_solve(
    struct d_simplex_workspace *ws,
    int entering_idx,
    double *solution)
{
    int k = ws->nreduced;
#if SIMPLEX_PROFILING
    acados_timer timer;
#endif

    for (int i = 0; i < k; ++i)
    {
        int row = ws->reduced_rows[i];
        BLASFEO_DVECEL(&ws->tmp_m_1, i) = entering_idx >= 0 ?
            d_simplex_A(ws, row, entering_idx) : ws->effective_rhs_data[row];
    }

    if (k > 0)
    {
#if SIMPLEX_PROFILING
        acados_tic(&timer);
#endif
        /* H x = b  <=>  T^T x = U^T L^T P x = b */
        blasfeo_dtrsv_utn(k, &ws->U_red, 0, 0, &ws->tmp_m_1, 0, &ws->tmp_m_1, 0);
        blasfeo_dtrsv_ltu(k, &ws->L_red, 0, 0, &ws->tmp_m_1, 0, &ws->tmp_m_1, 0);
        blasfeo_dvecin_sp(k, 1.0, &ws->tmp_m_1, 0, ws->row_perm, &ws->tmp_m_2, 0);
#if SIMPLEX_PROFILING
        ws->time_reduced_triangular_solves += acados_toc(&timer);
#endif

        for (int j = 0; j < k; ++j)
        {
            int basis_pos = ws->reduced_basis_pos[j];
            solution[basis_pos] = BLASFEO_DVECEL(&ws->tmp_m_2, j);
        }
    }

#if SIMPLEX_PROFILING
    acados_tic(&timer);
#endif
    for (int p = 0; p < ws->nsingleton; ++p)
    {
        int row = ws->singleton_rows[p];
        double value = entering_idx >= 0 ?
            d_simplex_A(ws, row, entering_idx) : ws->effective_rhs_data[row];
        BLASFEO_DVECEL(&ws->tmp_m_1, p) = value;
    }
    if (ws->nsingleton > 0 && k > 0)
        blasfeo_dgemv_n(ws->nsingleton, k, -1.0, &ws->G_mat, 0, 0,
                        &ws->tmp_m_2, 0, 1.0, &ws->tmp_m_1, 0,
                        &ws->tmp_m_1, 0);

    /* singleton columns are -e_row */
    for (int p = 0; p < ws->nsingleton; ++p)
        solution[ws->singleton_basis_pos[p]] = -BLASFEO_DVECEL(&ws->tmp_m_1, p);
#if SIMPLEX_PROFILING
    ws->time_reduced_gemv_recovery += acados_toc(&timer);
#endif
}

static void d_simplex_reduced_dual_solve(struct d_simplex_workspace *ws)
{
    int k = ws->nreduced;
    int ns = ws->nsingleton;
#if SIMPLEX_PROFILING
    acados_timer timer;
#endif

    for (int row = 0; row < ws->m; ++row)
        ws->lambda[row] = 0.0;

    /* lambda on the singleton rows is c_B / (-1); record its nonzeros */
    int nnz = 0, ns_nz = 0;
    for (int p = 0; p < ns; ++p)
    {
        int basis_pos = ws->singleton_basis_pos[p];
        int row = ws->singleton_rows[p];
        double value = -ws->cB[basis_pos]; /* singleton columns are -e_row */
        ws->lambda[row] = value;
        BLASFEO_DVECEL(&ws->tmp_m_2, p) = value;
        if (value != 0.0)
        {
            ws->lambda_nz[nnz++] = row;
            ws->sing_nz[ns_nz++] = p;
        }
    }

#if SIMPLEX_PROFILING
    acados_tic(&timer);
#endif
    for (int j = 0; j < k; ++j)
    {
        int basis_pos = ws->reduced_basis_pos[j];
        BLASFEO_DVECEL(&ws->tmp_m_1, j) = ws->cB[basis_pos];
    }
    if (ns_nz > 0 && k > 0)
    {
        if (4 * ns_nz < ns)
        {
            /* c_red -= G^T lambda_s over the nonzero singleton entries */
            for (int q = 0; q < ns_nz; ++q)
            {
                int p = ws->sing_nz[q];
                double lp = BLASFEO_DVECEL(&ws->tmp_m_2, p);
                for (int j = 0; j < k; ++j)
                    BLASFEO_DVECEL(&ws->tmp_m_1, j) -= BLASFEO_DMATEL(&ws->G_mat, p, j) * lp;
            }
        }
        else
        {
            blasfeo_dgemv_t(ns, k, -1.0, &ws->G_mat, 0, 0,
                            &ws->tmp_m_2, 0, 1.0, &ws->tmp_m_1, 0,
                            &ws->tmp_m_1, 0);
        }
    }
#if SIMPLEX_PROFILING
    ws->time_reduced_gemv_recovery += acados_toc(&timer);
#endif

    if (k > 0)
    {
#if SIMPLEX_PROFILING
        acados_tic(&timer);
#endif
        /* H^T y = c  <=>  T y = c  <=>  L U y = P c */
        blasfeo_dvecex_sp(k, 1.0, ws->row_perm, &ws->tmp_m_1, 0, &ws->tmp_m_2, 0);
        blasfeo_dtrsv_lnu(k, &ws->L_red, 0, 0, &ws->tmp_m_2, 0, &ws->tmp_m_2, 0);
        blasfeo_dtrsv_unn(k, &ws->U_red, 0, 0, &ws->tmp_m_2, 0, &ws->tmp_m_2, 0);
#if SIMPLEX_PROFILING
        ws->time_reduced_triangular_solves += acados_toc(&timer);
#endif

        for (int i = 0; i < k; ++i)
        {
            double value = BLASFEO_DVECEL(&ws->tmp_m_2, i);
            ws->lambda[ws->reduced_rows[i]] = value;
            if (value != 0.0)
                ws->lambda_nz[nnz++] = ws->reduced_rows[i];
        }
    }
    ws->lambda_nnz = nnz;
}

/*
 * Solve B^T lambda = c_B with the active factorization (reduced or full).
 *
 * Full basis: from P*B = L*U, B^T = U^T*L^T*P. First solve U^T*y = c_B and
 * L^T*z = y, then lambda = P^T*z.
 */
static void d_simplex_dual_solve(struct d_simplex_workspace *ws)
{
    int m = ws->m;
    if (ws->reduced_active)
    {
        d_simplex_reduced_dual_solve(ws);
        return;
    }
    ws->lambda_nnz = -1;

    blasfeo_pack_dvec(m, ws->cB, 1, &ws->tmp_m_2, 0);
    blasfeo_dtrsv_utn(m, &ws->LU_mat, 0, 0, &ws->tmp_m_2, 0, &ws->tmp_m_2, 0);
    blasfeo_dtrsv_ltu(m, &ws->LU_mat, 0, 0, &ws->tmp_m_2, 0, &ws->tmp_m_2, 0);

    /* lambda = P^T*z, i.e. lambda[row_perm[i]] = z[i]. */
    blasfeo_dvecin_sp(m, 1.0, &ws->tmp_m_2, 0, ws->row_perm, &ws->lu_work, 0);
    blasfeo_unpack_dvec(m, &ws->lu_work, 0, ws->lambda, 1);
}

static double d_simplex_variable_violation(
    const struct d_simplex_workspace *ws,
    int var,
    double value)
{
    double violation = 0.0;
    if (ws->has_lower[var])
        violation += fmax(0.0, ws->lower[var] - value);
    if (ws->has_upper[var])
        violation += fmax(0.0, value - ws->upper[var]);
    return violation;
}

/*
 * Phase-I cost of a basic variable: -1 below lower - tol, +1 above
 * upper + tol, 0 otherwise, with tol = feas_tol[var] (tol_infeas, times the
 * row norm for auxiliaries), so roundoff-level violations are not chased.
 */
static double d_simplex_phase1_cost(
    const struct d_simplex_workspace *ws,
    int var,
    double value)
{
    double tol = ws->feas_tol[var];
    if (ws->has_lower[var] && value < ws->lower[var] - tol)
        return -1.0;
    if (ws->has_upper[var] && value > ws->upper[var] + tol)
        return 1.0;
    return 0.0;
}

static int d_simplex_rebuild_effective_rhs(struct d_simplex_workspace *ws)
{
#if SIMPLEX_PROFILING
    acados_timer timer;
    acados_tic(&timer);
#endif
    const int n_nonbasic = ws->N - ws->m;
    blasfeo_dvecse(ws->nv, 0.0, &ws->rc, 0);
    for (int p = 0; p < n_nonbasic; ++p)
    {
        int var = ws->nonbasis_idx[p];
        if (var < ws->nv)
            blasfeo_dvecin1(ws->nonbasic_value[var], &ws->rc, var);
    }

    blasfeo_dgemv_n(ws->m, ws->nv, -1.0, &ws->Aphi, 0, 0,
                    &ws->rc, 0, 1.0, &ws->rhs, 0, &ws->tmp_m_1, 0);

    for (int p = 0; p < n_nonbasic; ++p)
    {
        int var = ws->nonbasis_idx[p];
        if (var < ws->nv)
            continue;
        int row = ws->singleton_row[var];
        BLASFEO_DVECEL(&ws->tmp_m_1, row) -=
            ws->singleton_value[var] * ws->nonbasic_value[var];
    }

    blasfeo_unpack_dvec(ws->m, &ws->tmp_m_1, 0, ws->effective_rhs_data, 1);
#if SIMPLEX_PROFILING
    ws->time_rhs_rebuild += acados_toc(&timer);
#endif
    ++ws->rhs_rebuild_count;
    /* x_B is recomputed from the fresh rhs, discarding incremental updates */
    ws->xB_valid = 0;
    return ACADOS_SUCCESS;
}

static void d_simplex_update_effective_rhs_column(
    struct d_simplex_workspace *ws,
    int var,
    double alpha)
{
    if (alpha == 0.0)
        return;
#if SIMPLEX_PROFILING
    acados_timer timer;
    acados_tic(&timer);
#endif

    int row = ws->singleton_row[var];
    if (row >= 0)
    {
        ws->effective_rhs_data[row] += alpha * ws->singleton_value[var];
    }
    else
    {
        for (int i = 0; i < ws->m; ++i)
            ws->effective_rhs_data[i] +=
                alpha * BLASFEO_DMATEL(&ws->Aphi, i, var);
    }
    ++ws->rhs_update_count;
#if SIMPLEX_PROFILING
    ws->time_rhs_update += acados_toc(&timer);
#endif
}

static int d_simplex_validate_event(
    const struct d_simplex_workspace *ws,
    int entering_idx,
    const double *x_B,
    const double *direction,
    const struct d_simplex_event *event)
{
    /* negated tests: also reject NaN */
    if (!(event->step >= 0.0) || !(event->entering_value == event->entering_value))
        return 0;

    const double tolerance = ws->tol_infeas *
        fmax(1.0, fabs(event->entering_value));
    if ((ws->has_lower[entering_idx] &&
         event->entering_value < ws->lower[entering_idx] - tolerance) ||
        (ws->has_upper[entering_idx] &&
         event->entering_value > ws->upper[entering_idx] + tolerance))
        return 0;

    if (event->kind == D_SIMPLEX_BOUND_FLIP)
    {
        /* the entering variable moves onto one of its bounds (the other one
         * for a variable at a bound, either one for a superbasic) */
        double expected_step = fabs(event->entering_value - ws->nonbasic_value[entering_idx]);
        int on_bound =
            (ws->has_lower[entering_idx] && event->entering_value == ws->lower[entering_idx]) ||
            (ws->has_upper[entering_idx] && event->entering_value == ws->upper[entering_idx]);
        return on_bound && fabs(event->step - expected_step) <= tolerance;
    }

    if (event->kind != D_SIMPLEX_BASIS_EXCHANGE)
        return 0;

    int var = ws->basis_idx[event->leaving_pos];
    double bound = event->leaving_side == 0 ?
        ws->lower[var] : ws->upper[var];
    double predicted = x_B[event->leaving_pos] -
        event->step * direction[event->leaving_pos];
    double scale = fmax(1.0, fmax(fabs(bound), fabs(predicted)));
    /* a leaving variable may sit anywhere inside its feasibility band */
    double leaving_tol = fmax(ws->tol_infeas * scale, ws->feas_tol[var]);
    if (!(fabs(direction[event->leaving_pos]) > ws->tol_pivot) ||
        !(fabs(predicted - bound) <= leaving_tol))
        return 0;

    return 1;
}

/*
 * Determines the pricing sign for nonbasic variable j with reduced cost rj.
 * Returns 1 if j is a pricing candidate, 0 if it must be skipped (fixed
 * variable), and -1 if the status is inconsistent (e.g. a basic variable
 * appearing in the nonbasic list).
 */
static int d_simplex_pricing_candidate_sign(
    const struct d_simplex_workspace *ws, int j, double rj, int *sign)
{
    if (ws->status[j] == D_SIMPLEX_AT_LOWER)
        *sign = 1;
    else if (ws->status[j] == D_SIMPLEX_AT_UPPER)
        *sign = -1;
    else if (ws->status[j] == D_SIMPLEX_FREE_NONBASIC ||
             ws->status[j] == D_SIMPLEX_SUPERBASIC)
        *sign = rj < 0.0 ? 1 : -1;
    else if (ws->status[j] == D_SIMPLEX_FIXED_NONBASIC)
        return 0;
    else
        return -1;
    return 1;
}

static struct d_simplex_event d_simplex_select_event(
    struct d_simplex_workspace *ws,
    int entering_idx,
    int entering_sign,
    double *x_B,
    double *direction,
    double directional_score)
{
    struct d_simplex_event event = {
        D_SIMPLEX_EVENT_FAILURE, -1, -1, 0.0, 0.0
    };
    const double direction_tol = ws->tol_pivot;
    const double pricing_tol = ws->tol;
    int event_count = 0;
    double slope = directional_score;
    double zero_step_correction = 0.0;
    int zero_step_leaving_pos = -1;
    int zero_step_leaving_side = -1;
    int zero_step_leaving_var = INT_MAX;

    /*
     * Right derivative of the banded Phase-I function along the move and the
     * breakpoints of the basic variables. cB holds the Phase-I cost at x_B,
     * so for a violated row the derivative cost * v equals the cB * v already
     * contained in directional_score: only feasible rows sitting in the band
     * of a bound and moving outward add to the slope. Bounds come from
     * lo_eff / up_eff (+-inf if absent: the breakpoint is then not finite).
     */
    const int m = ws->m;
    const int *basis_idx = ws->basis_idx;
    const double *lo_eff = ws->lo_eff, *up_eff = ws->up_eff, *band = ws->feas_tol;
    const double *cB = ws->cB;
    for (int i = 0; i < m; ++i)
    {
        int var = basis_idx[i];
        double value = x_B[i];
        double v = -direction[i];
        double lo = lo_eff[var], up = up_eff[var];
        int at_lower = fabs(value - lo) <= band[var];
        int at_upper = fabs(value - up) <= band[var];

        if (cB[i] == 0.0 && (at_lower || at_upper))
        {
            double correction = (at_lower ? fmax(0.0, -v) : 0.0) +
                                (at_upper ? fmax(0.0, v) : 0.0);
            slope += correction;
            zero_step_correction += correction;
            int active_side = at_lower && v < 0.0 ? 0 :
                (at_upper && v > 0.0 ? 1 : -1);
            if (active_side >= 0 && correction > 0.0 &&
                fabs(v) > direction_tol && var < zero_step_leaving_var)
            {
                zero_step_leaving_pos = i;
                zero_step_leaving_side = active_side;
                zero_step_leaving_var = var;
            }
        }

        if (fabs(v) <= direction_tol)
            continue;
        /* breakpoints of bounds inside the band are already accounted for
         * in the right derivative above */
        if (!at_lower)
        {
            double theta = (lo - value) / v;
            if (theta >= 0.0 && theta < INFINITY)
            {
                ws->thresholds[event_count].t = theta;
                ws->thresholds[event_count].alpha = fabs(v);
                ws->thresholds[event_count].row = i;
                ws->thresholds[event_count].side = 0;
                ++event_count;
            }
        }
        if (!at_upper)
        {
            double theta = (up - value) / v;
            if (theta >= 0.0 && theta < INFINITY)
            {
                ws->thresholds[event_count].t = theta;
                ws->thresholds[event_count].alpha = fabs(v);
                ws->thresholds[event_count].row = i;
                ws->thresholds[event_count].side = 1;
                ++event_count;
            }
        }
    }

    /* cap: distance of the entering variable to its own bound in the
     * direction of motion (other bound of a variable at a bound, either
     * bound of a superbasic one); reaching it first is a bound event */
    double cap = INFINITY;
    double cap_value = 0.0;
    int has_cap = 0;
    int entering_status = ws->status[entering_idx];
    if ((entering_status == D_SIMPLEX_AT_LOWER || entering_status == D_SIMPLEX_SUPERBASIC) &&
        entering_sign > 0 && ws->has_upper[entering_idx])
    {
        cap_value = ws->upper[entering_idx];
        cap = cap_value - ws->nonbasic_value[entering_idx];
        has_cap = 1;
    }
    else if ((entering_status == D_SIMPLEX_AT_UPPER || entering_status == D_SIMPLEX_SUPERBASIC) &&
             entering_sign < 0 && ws->has_lower[entering_idx])
    {
        cap_value = ws->lower[entering_idx];
        cap = ws->nonbasic_value[entering_idx] - cap_value;
        has_cap = 1;
    }
    if (has_cap && cap < 0.0)
        return event;


    if (slope >= -pricing_tol)
    {
        if (directional_score < -pricing_tol &&
            zero_step_correction > 0.0 && zero_step_leaving_pos >= 0)
        {
            event.kind = D_SIMPLEX_BASIS_EXCHANGE;
            event.leaving_pos = zero_step_leaving_pos;
            event.leaving_side = zero_step_leaving_side;
            event.entering_value = ws->nonbasic_value[entering_idx];
        }
        return event;
    }

    /*
     * Bound-flipping ratio test: visit the breakpoints in increasing order
     * until the slope of the Phase-I function becomes nonnegative. Usually
     * only a few breakpoints are visited, so they are kept in a binary heap
     * (O(n) to build, O(log n) per visited breakpoint) instead of sorting
     * all of them.
     */
    int heap_n = event_count;
    for (int i = heap_n / 2 - 1; i >= 0; --i)
        d_simplex_threshold_sift_down(ws->thresholds, heap_n, i);

    while (heap_n > 0)
    {
        double theta = ws->thresholds[0].t;
        double jump = 0.0;
        int selected_row = -1;
        int selected_side = -1;
        while (heap_n > 0 && ws->thresholds[0].t == theta)
        {
            struct d_simplex_threshold top;
            d_simplex_threshold_heap_pop(ws->thresholds, &heap_n, &top);
            jump += top.alpha;
            int candidate_row = top.row;
            int candidate_var = ws->basis_idx[candidate_row];
            int selected_var = selected_row >= 0 ?
                ws->basis_idx[selected_row] : INT_MAX;
            if (selected_row < 0 || candidate_var < selected_var ||
                (candidate_var == selected_var && top.side < selected_side))
            {
                selected_row = candidate_row;
                selected_side = top.side;
            }
        }

        if (theta == 0.0)
            continue;
        if (theta > 0.0 && slope < -pricing_tol && slope + jump >= -pricing_tol)
        {
            if (has_cap && cap < theta)
            {
                event.kind = D_SIMPLEX_BOUND_FLIP;
                event.step = cap;
                event.entering_value = cap_value;
                return event;
            }
            event.kind = D_SIMPLEX_BASIS_EXCHANGE;
            event.leaving_pos = selected_row;
            event.leaving_side = selected_side;
            event.step = theta;
            event.entering_value = ws->nonbasic_value[entering_idx] +
                entering_sign * theta;
            return event;
        }
        slope += jump;
    }

    if (has_cap)
    {
        event.kind = D_SIMPLEX_BOUND_FLIP;
        event.step = cap;
        event.entering_value = cap_value;
        return event;
    }
    return event;
}

/* Solve B x_B = effective rhs with the active factorization. */
static void d_simplex_compute_xB(struct d_simplex_workspace *ws)
{
    int m = ws->m;
    if (ws->reduced_active)
    {
        d_simplex_reduced_primal_solve(ws, -1, ws->x_B);
    }
    else
    {
        blasfeo_pack_dvec(m, ws->effective_rhs_data, 1, &ws->tmp_m_1, 0);
        blasfeo_dvecex_sp(m, 1.0, ws->row_perm, &ws->tmp_m_1, 0, &ws->lu_work, 0);
        blasfeo_dtrsv_lnu(m, &ws->LU_mat, 0, 0, &ws->lu_work, 0, &ws->lu_work, 0);
        blasfeo_dtrsv_unn(m, &ws->LU_mat, 0, 0, &ws->lu_work, 0, &ws->tmp_m_1, 0);
        blasfeo_unpack_dvec(m, &ws->tmp_m_1, 0, ws->x_B, 1);
    }
    ws->xB_solve_count++;
    ws->xB_valid = 1;
}

/* Solve B d = a_enter (column entering_idx of the Phase-I matrix) into ws->d_dir. */
static void d_simplex_direction_solve(struct d_simplex_workspace *ws, int entering_idx)
{
    int m = ws->m;
    if (ws->reduced_active)
    {
        d_simplex_reduced_primal_solve(ws, entering_idx, ws->d_dir);
    }
    else
    {
        d_simplex_column_extract(ws, entering_idx, &ws->tmp_m_1);
        blasfeo_dvecex_sp(m, 1.0, ws->row_perm, &ws->tmp_m_1, 0, &ws->lu_work, 0);
        blasfeo_dtrsv_lnu(m, &ws->LU_mat, 0, 0, &ws->lu_work, 0, &ws->lu_work, 0);
        blasfeo_dtrsv_unn(m, &ws->LU_mat, 0, 0, &ws->lu_work, 0, &ws->lu_work, 0);
        blasfeo_unpack_dvec(m, &ws->lu_work, 0, ws->d_dir, 1);
    }
}

/*
 * Commit a basis exchange: entering_idx (at nonbasis position entering_pos)
 * replaces the basic variable at event->leaving_pos, which leaves at the
 * bound event->leaving_side. Updates the effective rhs, the statuses, x_B
 * (x_B <- x_B - step * d_dir, with d_dir carrying the entering sign) and the
 * reduced representation / factors.
 */
static int d_simplex_commit_exchange(struct d_simplex_workspace *ws, int entering_idx,
                                     int entering_pos, const struct d_simplex_event *event,
                                     double tol_pivot)
{
    int m = ws->m;
    int leaving_pos = event->leaving_pos;
    int leaving_side = event->leaving_side;
    if (entering_pos < 0 || leaving_pos < 0 || leaving_pos >= m)
        return ACADOS_NAN_DETECTED;

    int leaving_idx = ws->basis_idx[leaving_pos];
    double leaving_bound = leaving_side == 0 ?
        ws->lower[leaving_idx] : ws->upper[leaving_idx];
    if (!(fabs(ws->d_dir[leaving_pos]) > tol_pivot))
        return ACADOS_NAN_DETECTED;

    int entering_is_structural = entering_idx < ws->nv;
    int leaving_is_structural = ws->singleton_row[leaving_idx] < 0;
    int pivot_type = entering_is_structural ?
        (leaving_is_structural ? 3 : 1) :
        (leaving_is_structural ? 2 : 0);
    ws->pivot_counts[pivot_type]++;

    // Commit the logical pivot and its RHS update before factor updates.
    double old_entering_value = ws->nonbasic_value[entering_idx];
    d_simplex_update_effective_rhs_column(
        ws, entering_idx, old_entering_value);
    d_simplex_update_effective_rhs_column(
        ws, leaving_idx, -leaving_bound);
    ws->nonbasic_value[entering_idx] = event->entering_value;
    ws->status[entering_idx] = D_SIMPLEX_BASIC;
    ws->nonbasic_value[leaving_idx] = leaving_bound;
    ws->status[leaving_idx] =
        ws->has_lower[leaving_idx] && ws->has_upper[leaving_idx] &&
        ws->lower[leaving_idx] == ws->upper[leaving_idx] ?
            D_SIMPLEX_FIXED_NONBASIC :
            (leaving_side == 0 ? D_SIMPLEX_AT_LOWER : D_SIMPLEX_AT_UPPER);
    ws->basis_idx[leaving_pos] = entering_idx;
    ws->nonbasis_idx[entering_pos] = leaving_idx;
    /* x_B <- x_B - step * d, and the entering variable takes the place
     * of the leaving one. The leaving variable is moved exactly onto its
     * bound (it may sit inside the feasibility band); that difference
     * is absorbed at the next rhs refresh. */
    for (int i = 0; i < m; ++i)
        ws->x_B[i] -= event->step * ws->d_dir[i];
    ws->x_B[leaving_pos] = event->entering_value;
    if (event->step == 0.0)
        ++ws->zero_step_count;
    ++ws->events_since_rhs_refresh;
    ws->rhs_fresh = 0;
    ws->numerical_retry_count = 0;

    int update_status = ACADOS_NAN_DETECTED;
    if (ws->reduced_rep_valid)
        update_status = d_simplex_update_reduced_representation(
            ws, pivot_type, entering_idx, leaving_idx, tol_pivot);

    ws->reduced_active = 0;
    if (update_status == ACADOS_SUCCESS)
    {
        ws->reduced_rep_valid = 1;
    }
    else
    {
        ws->reduced_rep_valid = 0;
        ws->reduced_lu_valid = 0;
    }
    return ACADOS_SUCCESS;
}

/*
 * Initialize the dual simplex phase on the slack start basis (all auxiliary
 * variables basic with zero cost, so lambda = 0 and d = c).
 */
static void d_simplex_dual_init(struct d_simplex_workspace *ws)
{
    unsigned int seed = 12345u;
    for (int j = 0; j < ws->N; ++j)
    {
        ws->dual_cost[j] = 0.0;
        if (j < ws->nv && (ws->status[j] == D_SIMPLEX_AT_LOWER ||
                           ws->status[j] == D_SIMPLEX_AT_UPPER))
        {
            double norm = ws->col_norm_sq[j];
            /* LCG constants of the ANSI C sample rand() implementation */
            seed = seed * 1103515245u + 12345u;
            double u = (double) ((seed >> 16) & 0x7fff) / 32768.0;
            double weight = (norm > 0.0 ? sqrt(norm) : 1.0) * (1.0 + 0.5 * u);
            ws->dual_cost[j] = ws->status[j] == D_SIMPLEX_AT_LOWER ? weight : -weight;
        }
        ws->dual_d[j] = ws->dual_cost[j];
    }
    for (int i = 0; i < ws->m; ++i)
        ws->dual_d[ws->basis_idx[i]] = 0.0;
}

/* Recompute the reduced costs d = c - Aphi^T lambda, B^T lambda = c_B,
 * to remove the drift of the incremental updates. Overwrites cB, lambda, rc. */
static int d_simplex_dual_refresh_reduced_costs(struct d_simplex_workspace *ws)
{
    for (int i = 0; i < ws->m; ++i)
        ws->cB[i] = ws->dual_cost[ws->basis_idx[i]];
    d_simplex_dual_solve(ws);
    d_simplex_compute_dual_product(ws);
    for (int k = 0; k < ws->N - ws->m; ++k)
    {
        int j = ws->nonbasis_idx[k];
        ws->dual_d[j] = ws->dual_cost[j] - d_simplex_column_dot_lambda(ws, j);
    }
    for (int i = 0; i < ws->m; ++i)
        ws->dual_d[ws->basis_idx[i]] = 0.0;
    return ACADOS_SUCCESS;
}

enum d_simplex_dual_result
{
    D_SIMPLEX_DUAL_PIVOTED,
    D_SIMPLEX_DUAL_INFEASIBLE,
    D_SIMPLEX_DUAL_FALLBACK,
    D_SIMPLEX_DUAL_NAN
};

static const char *d_simplex_dual_fallback_name(int reason)
{
    switch (reason)
    {
        case 0: return "none";
        case 1: return "iteration budget";
        case 2: return "certificate not confirmed";
        case 3: return "small or inconsistent pivot";
        case 4: return "solve failure";
        case 5: return "reduced-cost refresh failure";

    }
    return "unknown";
}

/*
 * One dual simplex iteration with leaving basic variable at position r
 * (primal infeasible beyond its band).
 *
 * Row: with c_r = -1 (x_r below lower) or +1 (above upper) and
 * lambda = B^{-T} c_r e_r, q = Aphi^T lambda. Increasing the violated side
 * of x_r requires moving nonbasic j with q_j > 0 up from its lower bound,
 * q_j < 0 down from its upper bound, or a free j in the direction of q_j.
 * Dual ratio test: the step t along the dual ray changes d_j -> d_j - t q_j;
 * candidate j becomes binding at t_j = |d_j| / |q_j|. The breakpoints are
 * visited in increasing order (heap); a boxed candidate whose full range
 * does not yet repair x_r is flipped to its other bound (dual BFRT). No
 * entering candidate means the dual is unbounded: lambda is a Farkas
 * certificate (sup over the boxes of q^T x < lambda^T rhs).
 */
static int d_simplex_dual_iteration(struct d_simplex_workspace *ws, int iter, int r,
                                    double tol_pivot)
{
    int m = ws->m;
    int N = ws->N;
    int leaving_var = ws->basis_idx[r];
    double x_r = ws->x_B[r];
    double c_r = d_simplex_phase1_cost(ws, leaving_var, x_r);
    if (c_r == 0.0)
        return D_SIMPLEX_DUAL_FALLBACK;
    int leaving_side = c_r < 0.0 ? 0 : 1;
    double leaving_bound = leaving_side == 0 ? ws->lower[leaving_var] : ws->upper[leaving_var];
    double s = -c_r; /* +1: x_r must increase, -1: decrease */

    /* row of B^{-1}: lambda = B^{-T} c_r e_r, q = Aphi^T lambda */
    for (int i = 0; i < m; ++i)
        ws->cB[i] = 0.0;
    ws->cB[r] = c_r;
    d_simplex_dual_solve(ws);
    d_simplex_compute_dual_product(ws);

    /* candidates and breakpoints */
    double q_max = 0.0;
    for (int k = 0; k < N - m; ++k)
        q_max = fmax(q_max, fabs(d_simplex_column_dot_lambda(ws, ws->nonbasis_idx[k])));
    double q_tol = SIMPLEX_DUAL_PIVOT_TOL * fmax(1.0, q_max);

    /* Free and superbasic variables have breakpoint t = 0. A free one with
     * a usable pivot (largest |q_j|) enters immediately; superbasic ones are
     * visited first in the bound-flipping test (key -|q_j|) and moved onto
     * their bound while x_r stays violated. */
    int count = 0;
    int zero_ratio_idx = -1, zero_ratio_pos = -1;
    double zero_ratio_q = 0.0;
    for (int k = 0; k < N - m; ++k)
    {
        int j = ws->nonbasis_idx[k];
        double q = d_simplex_column_dot_lambda(ws, j);
        if (fabs(q) <= q_tol)
            continue;
        double d = ws->dual_d[j];
        double ratio;
        if (ws->status[j] == D_SIMPLEX_AT_LOWER && q > 0.0)
            ratio = fmax(d, 0.0) / q;
        else if (ws->status[j] == D_SIMPLEX_AT_UPPER && q < 0.0)
            ratio = fmax(-d, 0.0) / -q;
        else if (ws->status[j] == D_SIMPLEX_SUPERBASIC)
            ratio = -fabs(q);
        else if (ws->status[j] == D_SIMPLEX_FREE_NONBASIC)
        {
            if (fabs(q) > zero_ratio_q)
            {
                zero_ratio_q = fabs(q);
                zero_ratio_idx = j;
                zero_ratio_pos = k;
            }
            continue;
        }
        else
            continue;
        ws->thresholds[count].t = ratio;
        ws->thresholds[count].alpha = fabs(q);
        ws->thresholds[count].row = j;
        ws->thresholds[count].side = k;
        ++count;
    }

    /* bound-flipping dual ratio test */
    double slope = s * (leaving_bound - x_r); /* remaining primal infeasibility of x_r */
    int heap_n = count;
    for (int i = heap_n / 2 - 1; i >= 0; --i)
        d_simplex_threshold_sift_down(ws->thresholds, heap_n, i);
    int entering_idx = -1, entering_pos = -1, flips = 0;
    double t_step = 0.0;
    if (zero_ratio_idx >= 0)
    {
        entering_idx = zero_ratio_idx;
        entering_pos = zero_ratio_pos;
        heap_n = 0;
    }
    while (heap_n > 0)
    {
        struct d_simplex_threshold top;
        d_simplex_threshold_heap_pop(ws->thresholds, &heap_n, &top);
        int j = top.row;
        t_step = fmax(top.t, 0.0);
        if (ws->status[j] == D_SIMPLEX_SUPERBASIC)
        {
            /* moving x_j in the direction sign(q_j) repairs x_r; it can move
             * by the distance to its bound on that side */
            double q = d_simplex_column_dot_lambda(ws, j);
            int up = q > 0.0;
            int bounded = up ? ws->has_upper[j] : ws->has_lower[j];
            double target = up ? ws->upper[j] : ws->lower[j];
            double range = bounded ? fabs(target - ws->nonbasic_value[j]) : INFINITY;
            if (bounded && slope - top.alpha * range > 0.0)
            {
                d_simplex_update_effective_rhs_column(ws, j, ws->nonbasic_value[j] - target);
                ws->status[j] = up ? D_SIMPLEX_AT_UPPER : D_SIMPLEX_AT_LOWER;
                ws->nonbasic_value[j] = target;
                slope -= top.alpha * range;
                ++flips;
                continue;
            }
            entering_idx = j;
            entering_pos = top.side;
            break;
        }
        int boxed = ws->has_lower[j] && ws->has_upper[j] &&
            ws->status[j] != D_SIMPLEX_FREE_NONBASIC;
        double range = boxed ? ws->upper[j] - ws->lower[j] : INFINITY;
        if (boxed && slope - top.alpha * range > 0.0)
        {
            /* flip j to its other bound: x_r moves by |q_j| * range */
            double old_value = ws->nonbasic_value[j];
            double new_value = ws->status[j] == D_SIMPLEX_AT_LOWER ? ws->upper[j] : ws->lower[j];
            d_simplex_update_effective_rhs_column(ws, j, old_value - new_value);
            ws->status[j] = ws->status[j] == D_SIMPLEX_AT_LOWER ?
                D_SIMPLEX_AT_UPPER : D_SIMPLEX_AT_LOWER;
            ws->nonbasic_value[j] = new_value;
            slope -= top.alpha * range;
            ++flips;
            continue;
        }
        entering_idx = j;
        entering_pos = top.side;
        break;
    }

    /* dual step: d_j <- d_j - t q_j on the nonbasic variables */
    for (int k = 0; k < N - m; ++k)
    {
        int j = ws->nonbasis_idx[k];
        ws->dual_d[j] -= t_step * d_simplex_column_dot_lambda(ws, j);
    }

    if (flips > 0)
    {
        ws->dual_flip_count += flips;
        ws->bound_flip_count += flips;
        ws->events_since_rhs_refresh += flips;
        ws->rhs_fresh = 0;
        ws->xB_valid = 0;
    }

    if (entering_idx < 0)
    {
        /* dual unbounded: lambda (with cB = c_r e_r) is a Farkas candidate */
        if (d_simplex_farkas_screen(ws, 1.0) &&
            d_simplex_try_farkas_certificate(ws, iter, 3))
            return D_SIMPLEX_DUAL_INFEASIBLE;
        ws->dual_fallback_reason = 2;
        return D_SIMPLEX_DUAL_FALLBACK;
    }

    if (flips > 0)
        d_simplex_compute_xB(ws);
    x_r = ws->x_B[r];
    if (d_simplex_phase1_cost(ws, leaving_var, x_r) != c_r)
    {
        /* the flips already repaired x_r (rounding): no basis change */
        return D_SIMPLEX_DUAL_PIVOTED;
    }

    /* primal step: x_r moves onto its violated bound */
    d_simplex_direction_solve(ws, entering_idx);
    double alpha_r = ws->d_dir[r];
    double theta = (x_r - leaving_bound) / alpha_r;
    int entering_sign = theta >= 0.0 ? 1 : -1;
    int expected_sign = ws->status[entering_idx] == D_SIMPLEX_AT_LOWER ? 1 :
        ws->status[entering_idx] == D_SIMPLEX_AT_UPPER ? -1 : entering_sign;
    if (!(fabs(alpha_r) > tol_pivot) || entering_sign != expected_sign)
    {
        ws->dual_fallback_reason = 3;
        return D_SIMPLEX_DUAL_FALLBACK;
    }
    for (int i = 0; i < m; ++i)
        ws->d_dir[i] *= entering_sign;

    struct d_simplex_event event;
    event.kind = D_SIMPLEX_BASIS_EXCHANGE;
    event.leaving_pos = r;
    event.leaving_side = leaving_side;
    event.step = fabs(theta);
    event.entering_value = ws->nonbasic_value[entering_idx] + theta;

    /* reduced costs of the exchanged pair */
    ws->dual_d[entering_idx] = 0.0;
    ws->dual_d[leaving_var] = s * t_step;

    if (d_simplex_commit_exchange(ws, entering_idx, entering_pos, &event, tol_pivot)
            != ACADOS_SUCCESS)
        return D_SIMPLEX_DUAL_NAN;
    ws->dual_iter_count++;
    return D_SIMPLEX_DUAL_PIVOTED;
}

/* Devex ranking measure of a pricing candidate (score < 0). */
static inline double d_simplex_pricing_measure(const struct d_simplex_workspace *ws,
                                               int j, double score)
{
#if SIMPLEX_DEVEX
    return score * score / ws->devex_w[j];
#else
    (void) ws; (void) j;
    return score * score;
#endif
}

#if SIMPLEX_DEVEX
static void d_simplex_devex_reset(struct d_simplex_workspace *ws)
{
    for (int j = 0; j < ws->N; ++j)
        ws->devex_w[j] = 1.0;
}

/*
 * Devex weight update for the basis exchange entering q at pivot row r,
 * before the factorization is updated (B is still the old basis). The pivot
 * row alpha_j = (B^{-1} a_j)_r = a_j^T rho, rho = B^{-T} e_r, is computed
 * with one dual solve and one (sparse) product. Overwrites cB, lambda, rc.
 */
static void d_simplex_devex_update(struct d_simplex_workspace *ws, int q, int r)
{
    const int m = ws->m, nv = ws->nv, N = ws->N;
    double alpha_q = ws->d_dir[r]; /* +-(B^{-1} a_q)_r, only its square is used */
    if (alpha_q == 0.0)
        return;
    for (int i = 0; i < m; ++i)
        ws->cB[i] = 0.0;
    ws->cB[r] = 1.0;
    d_simplex_dual_solve(ws);
    d_simplex_compute_dual_product(ws);

    double wq = ws->devex_w[q];
    double inv_alpha_q2 = 1.0 / (alpha_q * alpha_q);
    double max_w = 0.0;
    for (int k = 0; k < N - m; ++k)
    {
        int j = ws->nonbasis_idx[k];
        if (j == q)
            continue;
        double alpha_j = j < nv ? BLASFEO_DVECEL(&ws->rc, j) :
            ws->singleton_value[j] * ws->lambda[ws->singleton_row[j]];
        double candidate = alpha_j * alpha_j * inv_alpha_q2 * wq;
        if (candidate > ws->devex_w[j])
            ws->devex_w[j] = candidate;
        if (ws->devex_w[j] > max_w)
            max_w = ws->devex_w[j];
    }
    int leaving = ws->basis_idx[r];
    ws->devex_w[leaving] = fmax(wq * inv_alpha_q2, 1.0);
    if (max_w > SIMPLEX_DEVEX_RESET || ws->devex_w[leaving] > SIMPLEX_DEVEX_RESET)
    {
        d_simplex_devex_reset(ws);
#if SIMPLEX_PROFILING
        ++ws->devex_resets;
#endif
    }
}
#endif

int d_simplex_solve(
    struct d_dense_qp *qp,
    struct d_dense_qp_sol *qp_sol,
    struct d_simplex_workspace *ws)
{
    if (!ws || !ws->initialized)
        return ACADOS_QP_FAILURE;

    int m  = ws->m;
    int N  = ws->N;
    int nv = ws->nv;

    int max_iter = ws->max_iter > 0 ? ws->max_iter : 2000;

    double tol = ws->tol > 0.0 ? ws->tol : 1e-7;

    double tol_pivot = ws->tol_pivot > 0.0 ? ws->tol_pivot : 1e-12;
#if SIMPLEX_PROFILING
    acados_timer stage_timer;
#endif

#if SIMPLEX_PROFILING
    ws->time_rhs_rebuild = 0.0;
    ws->time_rhs_update = 0.0;
    ws->time_bound_flip = 0.0;
    ws->time_recovery_validation = 0.0;
    ws->time_factorization = 0.0;
    ws->time_basic_solve = 0.0;
    ws->time_phase1_objective = 0.0;
    ws->time_dual_solve = 0.0;
    ws->time_pricing = 0.0;
    ws->time_direction_solve = 0.0;
    ws->time_ratio_test = 0.0;
    ws->time_basis_update = 0.0;
    ws->time_reduced_maps = 0.0;
    ws->time_reduced_H_build = 0.0;
    ws->time_reduced_G_build = 0.0;
    ws->time_reduced_factorization_only = 0.0;
    ws->time_reduced_incremental_update = 0.0;
    ws->time_reduced_triangular_solves = 0.0;
    ws->time_reduced_gemv_recovery = 0.0;
    ws->time_farkas = 0.0;
#endif
    ws->rhs_rebuild_count = 0;
    ws->rhs_update_count = 0;
    ws->numerical_refresh_count = 0;
    ws->rhs_fresh = 1;
    ws->xB_valid = 0;
    ws->xB_solve_count = 0;
    ws->events_since_rhs_refresh = 0;
    ws->numerical_retry_count = 0;
    ws->reduced_H_scalar_reads = 0;
    ws->reduced_G_scalar_reads = 0;
    ws->lu_updates_since_refactor = 0;
    ws->lu_update_count = 0;
    ws->lu_update_fail_count = 0;
    ws->reduced_refactor_count = 0;
    ws->certificate_valid = 0;
    ws->certificate_margin = NAN;
    ws->certificate_scale = 0.0;
    ws->certificate_dual_residual = NAN;
    ws->certificate_dual_residual_scale = NAN;
    ws->certificate_check_count = 0;
    ws->farkas_dual_screen_count = 0;
    ws->farkas_row_screen_count = 0;
    ws->farkas_next_full_check_iter = 0;
    ws->farkas_detection = 0;
    ws->certificate_result = D_SIMPLEX_CERTIFICATE_NOT_FOUND;
    ws->certificate_diagnostic.failure = D_SIMPLEX_CERTIFICATE_FAILURE_NONE;
    ws->certificate_diagnostic.index = -1;
    ws->certificate_diagnostic.index_is_variable = 0;
    ws->certificate_diagnostic.value = NAN;
    ws->certificate_diagnostic.scale = NAN;
    ws->certificate_diagnostic.threshold = NAN;
    ws->reduced_k_min = INT_MAX;
    ws->reduced_k_max = 0;
    ws->reduced_k_samples = 0;
    ws->reduced_k_sum = 0.0;
    for (int i = 0; i < 4; ++i)
        ws->pivot_counts[i] = 0;
    ws->bound_flip_count = 0;
    ws->zero_step_count = 0;

    if (m == 0)
    {
        for (int i = 0; i < nv; i++)
            ws->w[i] = ws->nonbasic_value[i];

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif
        int solution_valid = d_simplex_validate_original_solution(ws);
        if (!solution_valid)
        {
#if SIMPLEX_PROFILING
            ws->time_recovery_validation += acados_toc(&stage_timer);
#endif
            return ACADOS_NAN_DETECTED;
        }

        blasfeo_pack_dvec(nv, ws->w, 1, qp_sol->v, 0);
        d_simplex_recover_soft_slacks(qp, qp_sol, ws);
#if SIMPLEX_PROFILING
        ws->time_recovery_validation += acados_toc(&stage_timer);
#endif

        ws->iter = 0;
        return ACADOS_SUCCESS;
    }

    if (ws->print_level > 1)
    {
        printf("\n================ SIMPLEX PHASE I START ================\n");
        printf("dims: nv=%d ne=%d ng=%d nb=%d m=%d N=%d "
               "nslack=%d\n",
               ws->nv, ws->ne, ws->ng, ws->nb,
               ws->m, ws->N, ws->nslack);

        print_idx_local(
            "initial basis_idx", ws->basis_idx, ws->m);
    }

    if (ws->print_level > 3)
    {
        printf("\nA = \n");
        blasfeo_print_dmat(ws->m, ws->nv, &ws->Aphi, 0, 0);

        printf("\nrhs = \n");
        blasfeo_print_dvec(ws->m, &ws->rhs, 0);

        printf("\nc = \n");
        blasfeo_print_dvec(ws->N, &ws->c, 0);
    }

    if (ws->print_level > 1)
        printf("======================================================\n");

    if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
        return ACADOS_NAN_DETECTED;

    ws->dual_mode = SIMPLEX_DUAL_DETECTION;
#if SIMPLEX_DEVEX
    d_simplex_devex_reset(ws);
#endif
#if SIMPLEX_PROFILING
    ws->devex_resets = 0;
    ws->time_devex = 0.0;
#endif
    ws->dual_iter_count = 0;
    ws->dual_step_count = 0;
    ws->dual_budget = -1;
    ws->dual_flip_count = 0;
    ws->dual_fallback_reason = 0;
    ws->dual_detected_infeasible = 0;
#if SIMPLEX_PROFILING
    ws->time_dual_simplex = 0.0;
#endif
    if (ws->dual_mode)
        d_simplex_dual_init(ws);

    for (int iter = 0; iter < max_iter; iter++)
    {
        ws->iter = iter;
        if (ws->events_since_rhs_refresh >= SIMPLEX_RHS_REFRESH_INTERVAL)
        {
            if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
                return ACADOS_NAN_DETECTED;
            ++ws->numerical_refresh_count;
            ws->events_since_rhs_refresh = 0;
            ws->rhs_fresh = 1;
        }

        if (ws->print_level > 1)
        {
            printf("\nSimplex iter %d:\n", iter);
            print_idx_local("\tbasis_idx", ws->basis_idx, m);
        }

        /*
         * STEP 1-2: Factorize the basis via LU, try first with reduced factors
         */

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        ws->reduced_active = 0;
        int reduced_status = ACADOS_QP_FAILURE;
        if (ws->reduced_rep_valid)
        {
            reduced_status = d_simplex_refresh_reduced_basis_positions(ws, tol_pivot);
            if (reduced_status == ACADOS_SUCCESS && !ws->reduced_lu_valid)
                reduced_status = d_simplex_refactor_reduced_H(ws, tol_pivot);
        }
        else
        {
            reduced_status = d_simplex_build_reduced_factorization(ws, tol_pivot);
        }

        if (reduced_status == ACADOS_SUCCESS)
        {
            ws->reduced_active = 1;
        }
        else
        {
            ws->reduced_rep_valid = 0;
            ws->reduced_lu_valid = 0;
            if (ws->print_level > 0)
                printf("d_simplex_solve: iter %d - full-basis refactorization.\n", iter);

            for (int j = 0; j < m; j++)
            {
                int bj = ws->basis_idx[j];
                if (bj < 0 || bj >= N)
                {
                    printf("d_simplex_solve: basis index %d out of bounds.\n",
                        bj);
                    return ACADOS_NAN_DETECTED;
                }
                d_simplex_column_extract(ws, bj, &ws->tmp_m_1);
                blasfeo_dcolin(m, &ws->tmp_m_1, 0, &ws->B_mat, 0, j);

            }

            blasfeo_dgetrf_rp(m, m, &ws->B_mat, 0, 0, &ws->LU_mat, 0, 0, ws->ipiv);

            /* BLASFEO ipiv is a sequence of row interchanges.
                Construct perm such that: (P*x)[i] = x[row_perm[i]]. */
            for (int i = 0; i < m; i++)
                ws->row_perm[i] = i;

            for (int i = 0; i < m; i++)
            {
                int pivot = ws->ipiv[i];

                if (pivot < 0 || pivot >= m)
                {
                    printf("d_simplex_solve: invalid LU pivot ipiv[%d]=%d.\n", i, pivot);
                    return ACADOS_NAN_DETECTED;
                }

                int tmp = ws->row_perm[i];
                ws->row_perm[i] = ws->row_perm[pivot];
                ws->row_perm[pivot] = tmp;
            }
            for (int i = 0; i < m; ++i)
            {
                if (!(fabs(BLASFEO_DMATEL(&ws->LU_mat, i, i)) > tol_pivot))
                    return ACADOS_NAN_DETECTED;
            }
        }

#if SIMPLEX_PROFILING
        ws->time_factorization += acados_toc(&stage_timer);
#endif

        if (ws->reduced_active)
        {
            int k = ws->nreduced;
            ws->reduced_k_sum += k;
            ws->reduced_k_samples++;
            if (k < ws->reduced_k_min)
                ws->reduced_k_min = k;
            if (k > ws->reduced_k_max)
                ws->reduced_k_max = k;
        }

        /* STEP 3: solve B x_B = rhs. */

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        /* x_B is updated incrementally after each pivot / bound flip; a
         * solve is only needed after the effective rhs was rebuilt (start,
         * periodic refresh, numerical retries, final feasibility check). */
        if (!ws->xB_valid)
            d_simplex_compute_xB(ws);

#if SIMPLEX_PROFILING
        ws->time_basic_solve += acados_toc(&stage_timer);
#endif

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        /* phi: sum of the violations beyond the feasibility band.
         * The point is feasible iff no basic variable lies outside the band,
         * i.e. max violation <= tol_infeas (normalized rows). */
        double phi = 0.0;
        double max_violation = 0.0;
        int max_violation_pos = -1;
        int num_violated = 0;
        for (int i = 0; i < m; i++)
        {
            int basis_var = ws->basis_idx[i];
            ws->cB[i] = d_simplex_phase1_cost(ws, basis_var, ws->x_B[i]);
            if (ws->cB[i] != 0.0)
            {
                ++num_violated;
                double violation = d_simplex_variable_violation(ws, basis_var, ws->x_B[i]);
                phi += violation;
                /* select by violation of the normalized constraint row */
                int vrow = ws->singleton_row[basis_var];
                if (vrow >= 0 && ws->row_norm[vrow] > 0.0)
                    violation /= ws->row_norm[vrow];
                if (violation > max_violation)
                {
                    max_violation = violation;
                    max_violation_pos = i;
                }
            }
        }

#if SIMPLEX_PROFILING
        ws->time_phase1_objective += acados_toc(&stage_timer);
#endif

        if (max_violation_pos < 0)
        {
            if (!ws->rhs_fresh)
            {
                if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
                    return ACADOS_NAN_DETECTED;
                ++ws->numerical_refresh_count;
                ws->rhs_fresh = 1;
                continue;
            }
#if SIMPLEX_PROFILING
            acados_tic(&stage_timer);
#endif

            for (int j = 0; j < nv; ++j)
                ws->w[j] = ws->nonbasic_value[j];
            for (int i = 0; i < m; ++i)
            {
                int var = ws->basis_idx[i];
                if (var < nv)
                    ws->w[var] = ws->x_B[i];
            }
            if (!d_simplex_validate_original_solution(ws))
            {
                ws->iter = iter;
#if SIMPLEX_PROFILING
                ws->time_recovery_validation += acados_toc(&stage_timer);
#endif
                return ACADOS_NAN_DETECTED;
            }
            blasfeo_pack_dvec(nv, ws->w, 1, qp_sol->v, 0);
            d_simplex_recover_soft_slacks(qp, qp_sol, ws);
#if SIMPLEX_PROFILING
            ws->time_recovery_validation += acados_toc(&stage_timer);
#endif
            ws->iter = iter;
            return ACADOS_SUCCESS;
        }

        if (ws->dual_mode)
        {
            if (ws->dual_budget < 0)
            {
                /* budget from the fraction of rows violated at the start */
                double violated_fraction = (double) num_violated / (m > 0 ? m : 1);
                ws->dual_budget = (int) ceil(SIMPLEX_DUAL_ITER_BUDGET * violated_fraction *
                                             violated_fraction * (m + nv));
                if (ws->dual_budget < SIMPLEX_DUAL_MIN_STEPS)
                    ws->dual_budget = SIMPLEX_DUAL_MIN_STEPS;
            }
            int dual_result = D_SIMPLEX_DUAL_FALLBACK;
#if SIMPLEX_PROFILING
            acados_tic(&stage_timer);
#endif
            if (ws->dual_step_count++ >= ws->dual_budget)
                ws->dual_fallback_reason = 1;
            else if (ws->dual_iter_count > 0 &&
                     ws->dual_iter_count % SIMPLEX_RHS_REFRESH_INTERVAL == 0 &&
                     d_simplex_dual_refresh_reduced_costs(ws) != ACADOS_SUCCESS)
                ws->dual_fallback_reason = 5;
            else
                dual_result = d_simplex_dual_iteration(ws, iter, max_violation_pos, tol_pivot);
#if SIMPLEX_PROFILING
            ws->time_dual_simplex += acados_toc(&stage_timer);
#endif
            if (ws->print_level > 1)
                printf("\tdual simplex iteration: result %d, phi=%.6e\n", dual_result, phi);
            if (dual_result == D_SIMPLEX_DUAL_INFEASIBLE)
            {
                ws->dual_detected_infeasible = 1;
                ws->iter = iter;
                return ACADOS_INFEASIBLE;
            }
            if (dual_result == D_SIMPLEX_DUAL_NAN)
                return ACADOS_NAN_DETECTED;
            if (dual_result == D_SIMPLEX_DUAL_PIVOTED)
                continue;
            /* fall back to primal Phase I from the current basis */
            ws->dual_mode = 0;
#if SIMPLEX_DEVEX
            d_simplex_devex_reset(ws);
#endif
            if (ws->print_level > 0)
                printf("d_simplex_solve: iter %d - dual simplex falls back to primal Phase I (%s).\n",
                       iter, d_simplex_dual_fallback_name(ws->dual_fallback_reason));
            if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
                return ACADOS_NAN_DETECTED;
            ++ws->numerical_refresh_count;
            ws->rhs_fresh = 1;
            continue;
        }

#if SIMPLEX_FARKAS_ROW_INTERVAL > 0
        /*
         * Early detection 2: single-row Farkas test on the most violated
         * basic variable x_r, with c_B = -+e_r and lambda = B^{-T} c_B. The
         * margin is l_r - max x_r (or min x_r - u_r) over the nonbasic
         * variables within their bounds. c_B is parked in d_dir, which is
         * unused before the direction solve.
         */
        if (iter % SIMPLEX_FARKAS_ROW_INTERVAL == 0)
        {
#if SIMPLEX_PROFILING
            acados_tic(&stage_timer);
#endif
            int r = max_violation_pos;
            for (int i = 0; i < m; ++i)
            {
                ws->d_dir[i] = ws->cB[i];
                ws->cB[i] = 0.0;
            }
            ws->cB[r] = ws->d_dir[r];

            d_simplex_dual_solve(ws);
            d_simplex_compute_dual_product(ws);
            ++ws->farkas_row_screen_count;
            int row_found = d_simplex_farkas_screen(ws, 1.0) &&
                d_simplex_try_farkas_certificate(ws, iter, 2);
            for (int i = 0; i < m; ++i)
                ws->cB[i] = ws->d_dir[i];
#if SIMPLEX_PROFILING
            ws->time_farkas += acados_toc(&stage_timer);
#endif
            if (row_found)
            {
                ws->iter = iter;
                return ACADOS_INFEASIBLE;
            }
        }
#endif

        /*
         * STEP 4: Solve B^T*lambda = c_B.
         */
#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        d_simplex_dual_solve(ws);

#if SIMPLEX_PROFILING
        ws->time_dual_solve += acados_toc(&stage_timer);
#endif

        if (ws->print_level > 1)
        {
            // check solution quality using inf norm of residuals
            for (int i = 0; i < m; i++)
            {
                double ri = -ws->effective_rhs_data[i];
                for (int j = 0; j < m; j++)
                {
                    ri += d_simplex_A(ws, i, ws->basis_idx[j]) * ws->x_B[j];
                }
                blasfeo_dvecin1(ri, &ws->tmp_m_1, i);
            }
            double res_Bx_inf;
            blasfeo_dvecnrm_inf(m, &ws->tmp_m_1, 0, &res_Bx_inf);
            printf("\tresidual: Bx - rhs = 0 ->  %.6e\n", res_Bx_inf);

            for (int i = 0; i < m; i++)
            {
                double ri = -ws->cB[i];
                for (int j = 0; j < m; j++)
                {
                    ri += d_simplex_A(ws, j, ws->basis_idx[i]) * ws->lambda[j];
                }
                blasfeo_dvecin1(ri, &ws->tmp_m_2, i);
            }

            double res_Btlam_inf;
            blasfeo_dvecnrm_inf(m, &ws->tmp_m_2, 0, &res_Btlam_inf);
            printf("\tresidual: B^{T} lambda - c_B = 0 ->  %.6e\n", res_Btlam_inf);

            // printf("\n\tVariable in the basis x_B = \n");
            // for (int i = 0; i < m; i++)
            //     // if (ws->basis_idx[i] < nv)
            //     printf("\t  %.6e", ws->x_B[i]);
            // printf("\n");
        }


        /*
         * STEP 5: block-partial pricing. The slack columns are -e_row, so
         * their reduced costs rj = lambda[singleton_row[j]] cost nothing and
         * all of them are scanned first. The best slack enters if its score
         * is below threshold = min(-0.1 phi, -1e-4). Otherwise, the
         * structural columns are priced in blocks of
         * SIMPLEX_PRICING_BLOCK_SIZE (one GEMV each), starting at block
         * iter mod n_blocks, until a block contains a score below -1e-4. The
         * best candidate seen enters. If m * nv <= SIMPLEX_FULL_PRICING_MAX_MN,
         * all columns are priced.
         */
        /* set if pricing left rc = -Aphi^T lambda for all structural columns */
        int rc_full = 0;

#if PARTIAL_PRICING

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        double threshold = fmin(-0.1 * phi, -1e-4);
        // printf("\tthreshold pricing: %.6e\n", threshold);

        /* Scores with magnitude below the dual tolerance opts->tol are
         * treated as zero: pivoting on them makes no meaningful progress, so
         * it is better to declare optimality and let the certificate check
         * decide between optimal and infeasible. */
        double pricing_zero_tol = tol;

        int entering_idx = -1;
        int entering_pos = -1;
        int entering_sign = 0;
        double most_neg_rc = -tol;
        double best_measure = 0.0;  /* pricing measure of the entering candidate */
        int slack_accepted = 0;

        /* 1. Slacks. If only slacks are nonbasic, the threshold is relaxed
         * to the scale of the infeasibility certificate. */
        double slack_best = -tol;
        double slack_measure = 0.0;
        int slack_idx = -1;
        int slack_pos = -1;
        int slack_sign = 0;
        int n_structural_nonbasic = 0;
        for (int k = 0; k < N - m; k++)
        {
            if (ws->nonbasis_idx[k] < nv)
                n_structural_nonbasic++;
        }
        int only_slacks_nonbasic = n_structural_nonbasic == 0;

        for (int k = 0; k < N - m; k++)
        {
            int j = ws->nonbasis_idx[k];
            if (j < nv)
                continue;

            double rj = ws->lambda[ws->singleton_row[j]];

            int candidate_sign;
            int sign_status = d_simplex_pricing_candidate_sign(ws, j, rj, &candidate_sign);
            if (sign_status < 0)
                return ACADOS_NAN_DETECTED;
            if (sign_status == 0)
                continue;

            double score = candidate_sign * rj;
            if (score > -pricing_zero_tol)
                continue;
            if (ws->print_level > 2)
                printf("\tPricing slack: j=%d sign=%d rc=%.6e score=%.6e\n",
                       j, candidate_sign, rj, score);
            double measure = d_simplex_pricing_measure(ws, j, score);
            if (measure > slack_measure ||
                (measure == slack_measure && j < slack_idx))
            {
                slack_measure = measure;
                slack_best = score;
                slack_idx = j;
                slack_pos = k;
                slack_sign = candidate_sign;
            }
        }

        /* small problems: full pricing, the best candidate overall wins */
        int full_pricing = (long) m * nv <= SIMPLEX_FULL_PRICING_MAX_MN;

        if (slack_idx >= 0)
        {
            double slack_threshold = only_slacks_nonbasic ?
                -fmax(ws->certificate_scale, 0.0) : threshold;
            /* The slack best is always kept as the provisional entering
             * candidate, even when it fails the threshold: structural
             * pricing may then only replace it with a strictly more
             * negative score, so the best overall candidate wins. */
            if (!full_pricing && slack_best < slack_threshold)
                slack_accepted = 1;
            entering_idx = slack_idx;
            entering_pos = slack_pos;
            entering_sign = slack_sign;
            most_neg_rc = slack_best;
            best_measure = slack_measure;
        }

        /* 2. Structural columns in blocks. The fixed threshold -1e-4 is used
         * here, since the dynamic one is not updated within an iteration and
         * could force pricing all blocks. */
        if (!slack_accepted)
        {
            blasfeo_pack_dvec(m, ws->lambda, 1, &ws->lu_work, 0);
            double struct_threshold = -1e-4;
            /* For small problems the whole GEMV is cheap: price all
             * structural variables at once instead of rotating blocks. */
            int n_blocks = full_pricing ? 1 :
                (nv + SIMPLEX_PRICING_BLOCK_SIZE - 1) / SIMPLEX_PRICING_BLOCK_SIZE;
            int first_block = full_pricing ? 0 : ws->iter % n_blocks;
            int blocks_priced = 0;
            for (int b = 0; b < n_blocks; b++)
            {
                ++blocks_priced;
                int block_idx = (first_block + b) % n_blocks;
                int col_offset = full_pricing ? 0 : block_idx * SIMPLEX_PRICING_BLOCK_SIZE;
                int block_size = full_pricing ? nv :
                    (nv - col_offset < SIMPLEX_PRICING_BLOCK_SIZE ?
                     nv - col_offset : SIMPLEX_PRICING_BLOCK_SIZE);

                d_simplex_lambda_times_A(ws, -1.0, col_offset, block_size, 1);

                double block_best = -tol;
                double block_measure = 0.0;
                int block_idx_entering = -1;
                int block_pos_entering = -1;
                int block_sign_entering = 0;

                for (int k = 0; k < N - m; k++)
                {
                    int j = ws->nonbasis_idx[k];
                    if (j < col_offset || j >= col_offset + block_size)
                        continue;

                    double rj = BLASFEO_DVECEL(&ws->rc, j);

                    int candidate_sign;
                    int sign_status = d_simplex_pricing_candidate_sign(ws, j, rj, &candidate_sign);
                    if (sign_status < 0)
                        return ACADOS_NAN_DETECTED;
                    if (sign_status == 0)
                        continue;

                    double score = candidate_sign * rj;
                    if (score > -pricing_zero_tol)
                        continue;
                    if (ws->print_level > 2)
                        printf("\tPricing structural: j=%d sign=%d rc=%.6e score=%.6e\n",
                               j, candidate_sign, rj, score);
                    double measure = d_simplex_pricing_measure(ws, j, score);
                    if (measure > block_measure ||
                        (measure == block_measure && j < block_idx_entering))
                    {
                        block_measure = measure;
                        block_best = score;
                        block_idx_entering = j;
                        block_pos_entering = k;
                        block_sign_entering = candidate_sign;
                    }
                }

                if (block_idx_entering >= 0 && block_measure > best_measure)
                {
                    entering_idx = block_idx_entering;
                    entering_pos = block_pos_entering;
                    entering_sign = block_sign_entering;
                    most_neg_rc = block_best;
                    best_measure = block_measure;
                    if (block_best < struct_threshold)
                        break;
                }
            }
            rc_full = blocks_priced == n_blocks;
        }

        if (ws->print_level > 1)
        {

            printf("\tentering_idx=%d (structural <%d) sign=%d directional_score=%.6e phi=%.6e\n",
                entering_idx, nv, entering_sign, most_neg_rc, phi);
            }

#if SIMPLEX_PROFILING
        ws->time_pricing += acados_toc(&stage_timer);
#endif

#else
        /*
         * STEP 5: Pricing.
         *
         * rc = c - Aphi^T*lambda.
         * lambda is packed into lu_work
         */

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        blasfeo_pack_dvec(m, ws->lambda, 1, &ws->lu_work, 0);
        // Structural columns are dense; every auxiliary column is -e_row.
        d_simplex_lambda_times_A(ws, -1.0, 0, nv, 1);
        rc_full = 1;

        int entering_idx = -1;
        int entering_pos = -1;
        int entering_sign = 0;
        double most_neg_rc = -tol;

        for (int k = 0; k < N - m; k++)
        {
            int j = ws->nonbasis_idx[k];

            if (j < 0 || j >= N)
            {
                printf("d_simplex_solve: invalid nonbasic index %d.\n",
                       j);
                return ACADOS_NAN_DETECTED;
            }

            int candidate_sign;
            if (ws->status[j] == D_SIMPLEX_AT_LOWER)
                candidate_sign = 1;
            else if (ws->status[j] == D_SIMPLEX_AT_UPPER)
                candidate_sign = -1;
            else if (ws->status[j] == D_SIMPLEX_FREE_NONBASIC ||
                     ws->status[j] == D_SIMPLEX_SUPERBASIC)
                candidate_sign = 0;
            else if (ws->status[j] == D_SIMPLEX_FIXED_NONBASIC)
                continue;
            else
                return ACADOS_NAN_DETECTED;

            double rj = j < nv ? BLASFEO_DVECEL(&ws->rc, j) :
                ws->lambda[ws->singleton_row[j]];

            if (candidate_sign == 0)
                candidate_sign = rj < 0.0 ? 1 : -1;

            double score = candidate_sign * rj;
            if (score > -tol)
                continue;
            if (ws->print_level > 2)
                printf("\tPricing: j=%d sign=%d rc=%.6e score=%.6e\n",
                       j, candidate_sign, rj, score);
            if (score < most_neg_rc ||
                (score == most_neg_rc && j < entering_idx))
            {
                most_neg_rc = score;
                entering_idx = j;
                entering_pos = k;
                entering_sign = candidate_sign;
            }
        }

        if (ws->print_level > 1)
            printf("\tentering_idx=%d sign=%d directional_score=%.6e\n",
                   entering_idx, entering_sign, most_neg_rc);

#if SIMPLEX_PROFILING
        ws->time_pricing += acados_toc(&stage_timer);
#endif
#endif

#if SIMPLEX_FARKAS_DUAL_INTERVAL > 0
        /*
         * Early detection 1: screen the Phase-I duals as a Farkas
         * certificate (weak-duality lower bound on phi). With full pricing,
         * rc = -Aphi^T lambda is already available and the screen is O(N);
         * otherwise the product is recomputed at a lower frequency.
         * When no entering variable exists, the terminal check below runs.
         */
        if (entering_idx >= 0 &&
            iter % (rc_full ? SIMPLEX_FARKAS_DUAL_INTERVAL :
                    SIMPLEX_FARKAS_DUAL_INTERVAL_PARTIAL) == 0)
        {
#if SIMPLEX_PROFILING
            acados_tic(&stage_timer);
#endif
            double q_sign = -1.0;
            if (!rc_full)
            {
                d_simplex_compute_dual_product(ws);
                q_sign = 1.0;
            }
            ++ws->farkas_dual_screen_count;
            int dual_found = d_simplex_farkas_screen(ws, q_sign) &&
                d_simplex_try_farkas_certificate(ws, iter, 1);
#if SIMPLEX_PROFILING
            ws->time_farkas += acados_toc(&stage_timer);
#endif
            if (dual_found)
            {
                ws->iter = iter;
                return ACADOS_INFEASIBLE;
            }
        }
#endif

        /*
         * STEP 6: Optimality/feasibility check.
         */
        if (entering_idx < 0)
        {
            d_simplex_compute_dual_residual(
                ws, &ws->certificate_dual_residual,
                &ws->certificate_dual_residual_scale);
            enum d_simplex_certificate_result certificate_result =
                d_simplex_run_certificate_check(ws);

            if (ws->print_level > 1)
            {
                printf("\tNo entering variable found.\n");
                printf("\tPhase-I objective phi = %.12e\n", phi);
                if (certificate_result != D_SIMPLEX_CERTIFICATE_VALID)
                    d_simplex_print_certificate_diagnostic(
                        certificate_result, &ws->certificate_diagnostic);
                printf("\nj=%d lower=%d %.17g upper=%d %.17g\n",
                    ws->certificate_diagnostic.index,
                    ws->has_lower[ws->certificate_diagnostic.index], ws->lower[ws->certificate_diagnostic.index],
                    ws->has_upper[ws->certificate_diagnostic.index], ws->upper[ws->certificate_diagnostic.index]);
            }

            if (ws->certificate_valid)
            {
                ws->iter = iter;
                return ACADOS_INFEASIBLE;
            }

            if (ws->numerical_retry_count++ < 2)
            {
                ws->reduced_rep_valid = 0;
                ws->reduced_lu_valid = 0;
                if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
                    return ACADOS_NAN_DETECTED;
                ++ws->numerical_refresh_count;
                ws->rhs_fresh = 1;
                continue;
            }
            /* Phase-I optimum with phi > 0 but no reliable certificate:
             * neither feasibility nor infeasibility could be established. */
            ws->iter = iter;
            return ACADOS_QP_FAILURE;
        }

        /*
        * STEP 7: Solve d = B^{-1}*a_enter.
        *
        * tmp_m_1 keeps the unpermuted entering column for the LU update.
        */
#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        d_simplex_direction_solve(ws, entering_idx);
        for (int i = 0; i < m; i++)
            ws->d_dir[i] *= entering_sign;

#if SIMPLEX_PROFILING
        ws->time_direction_solve += acados_toc(&stage_timer);
#endif

        /*
         * STEP 8: Ratio test.
         */
#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif

        struct d_simplex_event event = d_simplex_select_event(
            ws, entering_idx, entering_sign, ws->x_B, ws->d_dir,
            most_neg_rc);

#if SIMPLEX_PROFILING
        ws->time_ratio_test += acados_toc(&stage_timer);
#endif

        if (ws->print_level > 1)
        {
            printf("\tevent=%d leaving_pos=%d step=%.6e\n",
                   event.kind, event.leaving_pos, event.step);
        }

        if (event.kind == D_SIMPLEX_EVENT_FAILURE ||
            !d_simplex_validate_event(ws, entering_idx,
                ws->x_B, ws->d_dir, &event))
        {
            if (ws->numerical_retry_count++ < 2)
            {
                if (d_simplex_rebuild_effective_rhs(ws) != ACADOS_SUCCESS)
                    return ACADOS_NAN_DETECTED;
                ++ws->numerical_refresh_count;
                ws->rhs_fresh = 1;
                continue;
            }
            ws->iter = iter;
            return ACADOS_NAN_DETECTED;
        }

        if (event.kind == D_SIMPLEX_BOUND_FLIP)
        {
#if SIMPLEX_PROFILING
            acados_tic(&stage_timer);
#endif
            double old_value = ws->nonbasic_value[entering_idx];
            d_simplex_update_effective_rhs_column(
                ws, entering_idx, old_value - event.entering_value);
            /* the entering variable reached one of its bounds */
            ws->status[entering_idx] =
                ws->has_upper[entering_idx] &&
                event.entering_value == ws->upper[entering_idx] ?
                    D_SIMPLEX_AT_UPPER : D_SIMPLEX_AT_LOWER;
            ws->nonbasic_value[entering_idx] = event.entering_value;
            /* x_B <- x_B - step * d  (d_dir carries the entering sign) */
            for (int i = 0; i < m; ++i)
                ws->x_B[i] -= event.step * ws->d_dir[i];
            ++ws->bound_flip_count;
            ++ws->events_since_rhs_refresh;
            ws->rhs_fresh = 0;
            ws->numerical_retry_count = 0;
#if SIMPLEX_PROFILING
            ws->time_bound_flip += acados_toc(&stage_timer);
#endif
            continue;
        }

#if SIMPLEX_DEVEX
#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif
        d_simplex_devex_update(ws, entering_idx, event.leaving_pos);
#if SIMPLEX_PROFILING
        ws->time_devex += acados_toc(&stage_timer);
#endif
#endif

#if SIMPLEX_PROFILING
        acados_tic(&stage_timer);
#endif
        if (d_simplex_commit_exchange(ws, entering_idx, entering_pos, &event, tol_pivot)
                != ACADOS_SUCCESS)
            return ACADOS_NAN_DETECTED;
#if SIMPLEX_PROFILING
        ws->time_basis_update += acados_toc(&stage_timer);
#endif
    }

    /*
     * The iteration limit was reached while an entering variable was still
     * found. Run a final certificate check before reporting max-iter: the
     * throttled periodic checks may have missed a valid certificate.
     */
    {
        int infeasible_point = 0;
        for (int i = 0; i < m; i++)
        {
            int basis_var = ws->basis_idx[i];
            if (basis_var < 0 || basis_var >= N)
                return ACADOS_NAN_DETECTED;
            if (d_simplex_phase1_cost(ws, basis_var, ws->x_B[i]) != 0.0)
                infeasible_point = 1;
        }
        if (infeasible_point)
        {
            /* the certificate check needs the dual residual of the current
             * lambda; without it the result is always "uncertain" */
            d_simplex_compute_dual_residual(
                ws, &ws->certificate_dual_residual,
                &ws->certificate_dual_residual_scale);
            enum d_simplex_certificate_result certificate_result =
                d_simplex_run_certificate_check(ws);
            if (ws->print_level > 1)
                d_simplex_print_certificate_diagnostic(
                    certificate_result, &ws->certificate_diagnostic);
            if (certificate_result == D_SIMPLEX_CERTIFICATE_VALID)
            {
                ws->iter = max_iter;
                return ACADOS_INFEASIBLE;
            }
        }
    }

    ws->iter = max_iter;
    return ACADOS_MAXITER;
}

int dense_lp_feasibility(void *config, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    (void) config;
    (void) work_;

    dense_qp_in *qp_in = (dense_qp_in *) qp_in_;
    dense_qp_out *qp_out = (dense_qp_out *) qp_out_;
    dense_lp_feasibility_opts *opts = (dense_lp_feasibility_opts *) opts_;
    dense_lp_feasibility_memory *mem = (dense_lp_feasibility_memory *) mem_;

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, lp_timer;
    acados_tic(&tot_timer);

    struct d_simplex_workspace *sws = mem->simplex_workspace;
    int nv = qp_in->dim->nv;
    int ns = qp_in->dim->ns;
    blasfeo_dvecse(nv + 2 * ns, 0.0, qp_out->v, 0);

    sws->print_level = opts->print_level;
    sws->max_iter = opts->max_iter;
    sws->tol = opts->tol;
    sws->tol_infeas = opts->tol_infeas;
    sws->tol_pivot = opts->tol_pivot;

    acados_tic(&lp_timer);
    int st_setup = d_simplex_workspace_update_from_qp((struct d_dense_qp *) qp_in, sws);
    double lp_time_setup = acados_toc(&lp_timer);
    mem->time_lp_setup = lp_time_setup;

    if (st_setup != ACADOS_SUCCESS)
    {
        mem->status = st_setup;
        mem->iter = 0;
        mem->time_lp_solver_call = lp_time_setup;
        if (info)
        {
            info->solve_QP_time = lp_time_setup;
            info->interface_time = 0.0;
            info->total_time = acados_toc(&tot_timer);
            info->num_iter = 0;
            info->t_computed = 1;
        }

        return st_setup;
    }

    acados_tic(&lp_timer);
    int st_lp = d_simplex_solve((struct d_dense_qp *) qp_in, (struct d_dense_qp_sol *) qp_out, sws);
    double lp_time = acados_toc(&lp_timer);

    if (opts->print_level > 0)
    {
        printf("\n\n======== Simplex report ========");
        printf("\nVariables: %d original, %d original QP slacks.", sws->nv, qp_in->dim->ns);
        printf("\nPhase-I generic-constraint slacks (boxed): %d; equality residuals: %d.",
               sws->nineq_slack, sws->neq_residual);
        printf("\nRows: %d actual (%d equality, %d generic).",
               sws->m, sws->ne, sws->nineq_slack);
        printf("\nOriginal bound sets: %d.", sws->nb);
        printf("\nSetup time feas. LP: %f [sec]", lp_time_setup);
        printf("\nSolve time feas. LP: %f [sec]", lp_time);
        printf("\nTotal time feas. LP (setup + solve, reported as simplex_time): %f [sec]",
               lp_time_setup + lp_time);
#if SIMPLEX_PROFILING
        printf("\n  Effective-RHS rebuild: %f [sec] (%d rebuilds)",
               sws->time_rhs_rebuild, sws->rhs_rebuild_count);
        printf("\n  Effective-RHS updates: %f [sec] (%d column updates)",
               sws->time_rhs_update, sws->rhs_update_count);
        printf("\n  Bound flips: %f [sec] (%d flips)",
               sws->time_bound_flip, sws->bound_flip_count);
        printf("\n  Numerical refreshes: %d", sws->numerical_refresh_count);
        printf("\n  Basis factorization: %f [sec]", sws->time_factorization);
        printf("\n  Basic-variable solve: %f [sec] (%d solves)", sws->time_basic_solve,
               sws->xB_solve_count);
        printf("\n  Phase-I objective: %f [sec]", sws->time_phase1_objective);
        printf("\n  Dual-variable solve: %f [sec]", sws->time_dual_solve);
        printf("\n  Pricing: %f [sec]", sws->time_pricing);
        printf("\n  Direction solve: %f [sec]", sws->time_direction_solve);
        printf("\n  Ratio test: %f [sec]", sws->time_ratio_test);
        printf("\n  Basis update: %f [sec]", sws->time_basis_update);
        printf("\n  Reduced dimension k (avg/min/max): %.2f / %d / %d",
               sws->reduced_k_samples > 0 ? sws->reduced_k_sum / sws->reduced_k_samples : 0.0,
               sws->reduced_k_samples > 0 ? sws->reduced_k_min : 0,
               sws->reduced_k_max);
        printf("\n  Pivot count, singleton enters / singleton leaves: %d", sws->pivot_counts[0]);
        printf("\n  Pivot count, structural enters / singleton leaves: %d", sws->pivot_counts[1]);
        printf("\n  Pivot count, singleton enters / structural leaves: %d", sws->pivot_counts[2]);
        printf("\n  Pivot count, structural enters / structural leaves: %d", sws->pivot_counts[3]);
        printf("\n  Reduced map construction: %f [sec]", sws->time_reduced_maps);
        printf("\n  Reduced H construction: %f [sec]", sws->time_reduced_H_build);
        printf("\n  Reduced G construction: %f [sec]", sws->time_reduced_G_build);
        printf("\n  Reduced factorization only: %f [sec]", sws->time_reduced_factorization_only);
        printf("\n  Reduced incremental update: %f [sec]", sws->time_reduced_incremental_update);
        printf("\n  Reduced LU: %d refactorizations, %d updates, %d failed updates",
               sws->reduced_refactor_count, sws->lu_update_count, sws->lu_update_fail_count);
        printf("\n  Reduced H/G scalar reads: %d / %d",
             sws->reduced_H_scalar_reads, sws->reduced_G_scalar_reads);
        printf("\n  Reduced triangular solves: %f [sec]", sws->time_reduced_triangular_solves);
        printf("\n  Reduced GEMV/recovery: %f [sec]", sws->time_reduced_gemv_recovery);
        printf("\n  Dual simplex: %f [sec] (%d pivots, %d bound flips, budget %d steps, fallback: %s)",
               sws->time_dual_simplex, sws->dual_iter_count, sws->dual_flip_count,
               sws->dual_budget, d_simplex_dual_fallback_name(sws->dual_fallback_reason));
        printf("\n  Devex weight updates: %f [sec] (%d resets)", sws->time_devex, sws->devex_resets);
        printf("\n  Early Farkas screens: %f [sec] (%d dual, %d single-row)",
               sws->time_farkas, sws->farkas_dual_screen_count,
               sws->farkas_row_screen_count);
        printf("\n  Recovery and validation: %f [sec]", sws->time_recovery_validation);
        double profiled_time = sws->time_factorization + sws->time_basic_solve +
            sws->time_phase1_objective + sws->time_dual_solve + sws->time_pricing +
            sws->time_direction_solve + sws->time_ratio_test +
            sws->time_basis_update + sws->time_bound_flip +
            sws->time_rhs_rebuild + sws->time_recovery_validation +
            sws->time_farkas + sws->time_dual_simplex + sws->time_devex;
        printf("\n  Other simplex work: %f [sec]", fmax(0.0, lp_time - profiled_time));
#endif
        printf("\nSimplex status: %d", st_lp);
        printf("\nSimplex iter: %d", sws->iter);
        printf("\nCertificate tolerance: %.3e (relative)", SIMPLEX_CERTIFICATE_TOL);
        printf("\nDual residual infinity norm: %.12e", sws->certificate_dual_residual);
        printf("\nDual residual scale: %.12e", sws->certificate_dual_residual_scale);
        printf("\nRelative dual residual: %.12e",
               sws->certificate_dual_residual_scale > 0.0 ?
                   sws->certificate_dual_residual /
                       sws->certificate_dual_residual_scale :
                   sws->certificate_dual_residual);
        const char *certificate_status =
            sws->certificate_result == D_SIMPLEX_CERTIFICATE_VALID ?
                "valid" :
            sws->certificate_result == D_SIMPLEX_CERTIFICATE_NUMERICALLY_UNCERTAIN ?
                "numerically uncertain" : "not found";
        printf("\nInfeasibility certificate: %s", certificate_status);
        if (sws->certificate_valid)
            printf(" (found by %s)",
                   sws->farkas_detection == 1 ? "early Phase-I dual screen" :
                   sws->farkas_detection == 2 ? "early single-row test" :
                   sws->farkas_detection == 3 ? "dual simplex (unbounded dual ray)" :
                   "terminal Phase-I check");
        if (sws->certificate_valid || opts->print_level > 1)
        {
            printf("\nCertificate margin: %.12e", sws->certificate_margin);
            printf("\nCertificate scale: %.12e", sws->certificate_scale);
            printf("\nCertificate checks: %d", sws->certificate_check_count);
            printf("\nCertificate failure reason: %s",
                   d_simplex_certificate_failure_name(
                       sws->certificate_diagnostic.failure));
            if (sws->certificate_diagnostic.index >= 0)
            {
                printf("\nCertificate diagnostic %s: %d",
                       sws->certificate_diagnostic.index_is_variable ?
                           "variable" : "row",
                       sws->certificate_diagnostic.index);
                if (sws->certificate_diagnostic.index_is_variable)
                    printf("\nCertificate q: %.17e"
                           "\nCertificate q cancellation scale: %.17e"
                           "\nCertificate q sign tolerance: %.17e",
                           sws->certificate_diagnostic.value,
                           sws->certificate_diagnostic.scale,
                           sws->certificate_diagnostic.threshold);
            }
            else if (sws->certificate_diagnostic.failure !=
                     D_SIMPLEX_CERTIFICATE_FAILURE_NONE)
            {
                printf("\nCertificate diagnostic value: %.17e"
                       "\nCertificate diagnostic scale: %.17e"
                       "\nCertificate diagnostic threshold: %.17e",
                       sws->certificate_diagnostic.value,
                       sws->certificate_diagnostic.scale,
                       sws->certificate_diagnostic.threshold);
            }
        }
        printf("\nBound flips: %d; zero-step pivots: %d.",
               sws->bound_flip_count, sws->zero_step_count);
        printf("\nLP primal solution: \n");
        blasfeo_print_tran_dvec(qp_in->dim->nv + 2 * qp_in->dim->ns, qp_out->v, 0);
        printf("===============================\n");

    }

    // ------------------------------------------------------------
    // 3) fill bookkeeping
    // ------------------------------------------------------------
    /* the setup (building the Phase-I LP from the dense QP) is part of the
     * method, so it is included in the reported LP time */
    mem->time_lp_solver_call = lp_time_setup + lp_time;
    mem->iter = sws->iter;  // if available in your ws, otherwise set 0
    mem->status = st_lp;
    if (info)
    {
        info->solve_QP_time = lp_time;     // LP feasibility solve time
        info->interface_time = lp_time_setup; // building the Phase-I LP
        info->total_time = acados_toc(&tot_timer);
        info->num_iter = mem->iter;
        info->t_computed = 1;
    }

    return st_lp;
}