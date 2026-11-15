/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// blasfeo
#include "blasfeo_d_aux.h"

// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_daqp.h"
#include "acados/dense_qp/dense_qp_daqp.h"
#include "acados/utils/math.h"
#include "acados/utils/mem.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// daqp
#include "daqp/include/api.h"
#include "daqp/include/daqp.h"
#include "daqp/include/auxiliary.h"
#include "daqp/include/multistage.h"



/************************************************
 * opts
 ************************************************/

acados_size_t ocp_qp_daqp_opts_calculate_size(void *config_, void *dims_)
{
    acados_size_t size = 0;
    size += sizeof(ocp_qp_daqp_opts);
    size += sizeof(DAQPSettings);
    size += 8;
    return size;
}



void *ocp_qp_daqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    ocp_qp_daqp_opts *opts = (ocp_qp_daqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_qp_daqp_opts);

    align_char_to(8, &c_ptr);
    opts->daqp_opts = (DAQPSettings *) c_ptr;
    c_ptr += sizeof(DAQPSettings);

    assert((char *) raw_memory + ocp_qp_daqp_opts_calculate_size(config_, dims_) >= c_ptr);

    return (void *) opts;
}



void ocp_qp_daqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_qp_daqp_opts *opts = opts_;
    daqp_default_settings(opts->daqp_opts);
    opts->warm_start = 1;
    opts->print_level = 0;
}



void ocp_qp_daqp_opts_update(void *config_, void *dims_, void *opts_)
{
    return;
}



void ocp_qp_daqp_opts_set(void *config_, void *opts_, const char *field, void *value)
{
    ocp_qp_daqp_opts *opts = opts_;
    if (!strcmp(field, "tol_stat") || !strcmp(field, "tol_eq") || !strcmp(field, "tol_comp"))
    {
        // stationarity, equality constraints, and complementarity are
        // satisfied by the working set
    }
    else if (!strcmp(field, "tol_ineq"))
    {
        // At most the default of DAQP: a looser tolerance gains little for an
        // active-set method, but accepts violations of soft constraints that
        // the L1 penalty of an NLP solver weighs heavily
        double *tol = value;
        opts->daqp_opts->primal_tol = MIN(*tol, DAQP_DEFAULT_PRIM_TOL);
    }
    else if (!strcmp(field, "iter_max"))
    {
        int *iter_max = value;
        opts->daqp_opts->iter_limit = *iter_max;
    }
    else if (!strcmp(field, "warm_start"))
    {
        int *warm_start = value;
        opts->warm_start = *warm_start;
    }
    else if (!strcmp(field, "print_level"))
    {
        int *print_level = value;
        opts->print_level = *print_level;
    }
    else
    {
        printf("\nerror: ocp_qp_daqp_opts_set: wrong field: %s\n", field);
        exit(1);
    }
}



void ocp_qp_daqp_opts_get(void *config_, void *opts_, const char *field, void *value)
{
    printf("\nerror: ocp_qp_daqp_opts_get: not implemented for field %s\n", field);
    exit(1);
}



/************************************************
 * memory
 ************************************************/

// Number of arrays of pointers (to stage data) in DAQPMultistage that are used
#define OCP_QP_DAQP_NPTR 7

// Sizes of the data in the format of DAQP
static void ocp_qp_daqp_sizes(const ocp_qp_dims *dims, acados_size_t *nd, acados_size_t *ni,
        int *nz, int *m, int *npi)
{
    int N = dims->N;
    int *nx = dims->nx, *nu = dims->nu, *nb = dims->nb, *ng = dims->ng;
    *nd = 0; *ni = 0; *nz = 0; *m = 0; *npi = 0;
    for (int t = 0; t <= N; t++)
    {
        int nx1 = t < N ? nx[t+1] : 0, nut = t < N ? nu[t] : 0;
        *nd += nx1 + nx[t] + nut;  // b, q, r
        *nd += 2*nb[t] + 2*ng[t];  // lb, ub, lg, ug
        *ni += nb[t];  // idxb
        *nz += nx[t] + nut;
        *m += nb[t] + ng[t];
        *npi += nx1;
    }
}

// Size of the copies of the matrices (the dimensions of qp_in, so that the
// panels have the same layout)
static acados_size_t ocp_qp_daqp_copies_size(const ocp_qp_dims *dims)
{
    int N = dims->N;
    acados_size_t size = 3*(N+1)*sizeof(struct blasfeo_dmat);
    for (int t = 0; t <= N; t++)
    {
        int nux = (t < N ? dims->nu[t] : 0) + dims->nx[t];
        if (t < N)
            size += blasfeo_memsize_dmat(nux+1, dims->nx[t+1]);
        size += blasfeo_memsize_dmat(nux+1, nux);
        size += blasfeo_memsize_dmat(nux, dims->ng[t]);
    }
    return size + 64;
}


// Size of the memory of DAQP's multistage solver (from the dimensions alone), with the number of variables n and soft
// constraints ns of its core workspace
static acados_size_t ocp_qp_daqp_block_size(const ocp_qp_dims *dims, int *n, int *ns)
{
    DAQPMultistage ocp;
    memset(&ocp, 0, sizeof(ocp));
    ocp.N = dims->N; ocp.nx = dims->nx; ocp.nu = dims->nu; ocp.nb = dims->nb; ocp.ng = dims->ng;
    *ns = 0;  // (each slack softens a single constraint)
    for (int t = 0; t <= dims->N; t++)
        *ns += dims->ns[t];
    *n = daqp_multistage_nvmax(&ocp);
    return daqp_multistage_memsize(&ocp, *ns);
}



acados_size_t ocp_qp_daqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_qp_dims *dims = dims_;
    int N = dims->N, nz, m, npi, n, ns;
    acados_size_t nd, ni;
    ocp_qp_daqp_sizes(dims, &nd, &ni, &nz, &m, &npi);

    acados_size_t size = 0;
    size += sizeof(ocp_qp_daqp_memory);
    size += ocp_qp_daqp_block_size(dims, &n, &ns);
    size += daqp_workspace_calculate_size(n, m, 0, ns, 0);
    size += (OCP_QP_DAQP_NPTR*(N+1)) * sizeof(double *);
    size += (N+1) * sizeof(int *);  // idxb
    size += (nd + nz + m + npi) * sizeof(double);
    size += ocp_qp_daqp_copies_size(dims);
    size += 6 * m * sizeof(double);  // rho_l, rho_u, w_l, w_u, s_l, s_u
    size += (ni + 4*(N+1) + m) * sizeof(int);  // idxb, nx, nu, nb, ng, sense
    size += 2*8;
    return size;
}



void *ocp_qp_daqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_qp_dims *dims = dims_;
    int N = dims->N, nz, m, npi, n, ns;
    int *nx = dims->nx, *nu = dims->nu, *nb = dims->nb, *ng = dims->ng;
    acados_size_t nd, ni;
    ocp_qp_daqp_sizes(dims, &nd, &ni, &nz, &m, &npi);

    char *c_ptr = (char *) raw_memory;

    ocp_qp_daqp_memory *mem = (ocp_qp_daqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_qp_daqp_memory);
    memset(mem, 0, sizeof(*mem));

    // DAQP (nothing is allocated by its multistage solver)
    mem->daqp_mem_size = ocp_qp_daqp_block_size(dims, &n, &ns);
    mem->daqp_mem = c_ptr;
    c_ptr += mem->daqp_mem_size;
    align_char_to(8, &c_ptr);
    mem->work = daqp_workspace_assign(n, m, 0, ns, 0, c_ptr);
    c_ptr += daqp_workspace_calculate_size(n, m, 0, ns, 0);

    // arrays of pointers
    DAQPMultistage *ocp = &mem->ocp;
    memset(ocp, 0, sizeof(*ocp));
    ocp->N = N;
    // (the matrices are BLASFEO matrices of qp_in: A, B, Q, S, R, C, D stay NULL)
    double ***fam[OCP_QP_DAQP_NPTR] = {&ocp->b, &ocp->q, &ocp->r, &ocp->lb, &ocp->ub, &ocp->lg, &ocp->ug};
    for (int k = 0; k < OCP_QP_DAQP_NPTR; k++)
    {
        *fam[k] = (double **) c_ptr;
        c_ptr += (N+1) * sizeof(double *);
    }
    ocp->idxb = (int **) c_ptr;
    c_ptr += (N+1) * sizeof(int *);

    // doubles
    for (int t = 0; t <= N; t++)
    {
        int nx1 = t < N ? nx[t+1] : 0, nut = t < N ? nu[t] : 0;
        ocp->b[t] = (double *) c_ptr; c_ptr += nx1 * sizeof(double);
        ocp->q[t] = (double *) c_ptr; c_ptr += nx[t] * sizeof(double);
        ocp->r[t] = (double *) c_ptr; c_ptr += nut * sizeof(double);
        ocp->lb[t] = (double *) c_ptr; c_ptr += nb[t] * sizeof(double);
        ocp->ub[t] = (double *) c_ptr; c_ptr += nb[t] * sizeof(double);
        ocp->lg[t] = (double *) c_ptr; c_ptr += ng[t] * sizeof(double);
        ocp->ug[t] = (double *) c_ptr; c_ptr += ng[t] * sizeof(double);
    }
    mem->z = (double *) c_ptr; c_ptr += nz * sizeof(double);
    mem->lam = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->pi = (double *) c_ptr; c_ptr += npi * sizeof(double);
    mem->rho_l = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->rho_u = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->w_l = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->w_u = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->s_l = (double *) c_ptr; c_ptr += m * sizeof(double);
    mem->s_u = (double *) c_ptr; c_ptr += m * sizeof(double);

    // ints
    ocp->nx = (int *) c_ptr; c_ptr += (N+1) * sizeof(int);
    ocp->nu = (int *) c_ptr; c_ptr += (N+1) * sizeof(int);
    ocp->nb = (int *) c_ptr; c_ptr += (N+1) * sizeof(int);
    ocp->ng = (int *) c_ptr; c_ptr += (N+1) * sizeof(int);
    mem->sense = (int *) c_ptr; c_ptr += m * sizeof(int);
    ocp->sense = mem->sense;
    for (int i = 0; i < m; i++)
    {
        mem->sense[i] = 0;
        mem->rho_l[i] = mem->rho_u[i] = mem->w_l[i] = mem->w_u[i] = 0;
        mem->s_l[i] = mem->s_u[i] = 0;
    }
    for (int t = 0; t <= N; t++)
    {
        ocp->nx[t] = nx[t];
        ocp->nu[t] = t < N ? nu[t] : 0;
        ocp->nb[t] = nb[t];
        ocp->ng[t] = ng[t];
        ocp->idxb[t] = (int *) c_ptr; c_ptr += nb[t] * sizeof(int);
    }

    // copies of the matrices
    align_char_to(8, &c_ptr);
    mem->cBAbt = (struct blasfeo_dmat *) c_ptr;
    c_ptr += 3*(N+1)*sizeof(struct blasfeo_dmat);
    mem->cRSQrq = mem->cBAbt + N+1;
    mem->cDCt = mem->cRSQrq + N+1;
    align_char_to(64, &c_ptr);
    for (int t = 0; t <= N; t++)
    {
        int nux = (t < N ? nu[t] : 0) + nx[t];
        if (t < N)
            assign_and_advance_blasfeo_dmat_mem(nux+1, nx[t+1], mem->cBAbt+t, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nux+1, nux, mem->cRSQrq+t, &c_ptr);
        assign_and_advance_blasfeo_dmat_mem(nux, ng[t], mem->cDCt+t, &c_ptr);
    }

    // (the data is compared with the previous data at each call)
    memset(ocp->b[0], 0, nd * sizeof(double));
    mem->setup_done = 0;

    assert((char *) raw_memory + ocp_qp_daqp_memory_calculate_size(config_, dims, opts_) >= c_ptr);

    return mem;
}



void ocp_qp_daqp_memory_get(void *config_, void *mem_, const char *field, void* value)
{
    ocp_qp_daqp_memory *mem = mem_;

    if (!strcmp(field, "time_qp_solver_call"))
    {
        double *tmp_ptr = value;
        *tmp_ptr = mem->time_qp_solver_call;
    }
    else if (!strcmp(field, "iter"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->iter;
    }
    else if (!strcmp(field, "status"))
    {
        int *tmp_ptr = value;
        *tmp_ptr = mem->status;
    }
    else
    {
        printf("\nerror: ocp_qp_daqp_memory_get: field %s not available\n", field);
        exit(1);
    }
}



// Drop the setup of DAQP's multistage solver (its memory is in mem)
static void ocp_qp_daqp_free_work(ocp_qp_daqp_memory *mem)
{
    if (mem->setup_done)
    {
        free_daqp_multistage(mem->work);
        mem->setup_done = 0;
    }
}



void ocp_qp_daqp_memory_reset(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    // a cold start at the next call
    ocp_qp_daqp_free_work(mem_);
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_qp_daqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    return 0;
}



/************************************************
 * functions
 ************************************************/

// sC <-- sA (the first m rows and n columns). Returns 1 if it has changed
static int ocp_qp_daqp_sync_mat(int m, int n, struct blasfeo_dmat *sA, struct blasfeo_dmat *sC)
{
    int changed = 0;
    if (m <= 0 || n <= 0)
        return 0;
#if defined(MF_PANELMAJ)
    if (sA->cn == sC->cn)
    {
        // panels of D_PS rows (of cn columns each, contiguous): the full panels
        // at once, and the first n columns of the last one
        const int ps = D_PS, cn = sA->cn, mf = m/ps*ps, rows = m - mf;
        double *a = sA->pA + mf*cn, *c = sC->pA + mf*cn;
        if (mf > 0 && memcmp(sA->pA, sC->pA, mf*cn*sizeof(double)) != 0)
        {
            memcpy(sC->pA, sA->pA, mf*cn*sizeof(double));
            changed = 1;
        }
        for (int j = 0; rows > 0 && j < n; j++)
            if (memcmp(a+j*ps, c+j*ps, rows*sizeof(double)) != 0)
            {
                memcpy(c+j*ps, a+j*ps, rows*sizeof(double));
                changed = 1;
            }
        return changed;
    }
#endif
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            if (BLASFEO_DMATEL(sA, i, j) != BLASFEO_DMATEL(sC, i, j))
            {
                BLASFEO_DMATEL(sC, i, j) = BLASFEO_DMATEL(sA, i, j);
                changed = 1;
            }
    return changed;
}

// v = n elements of sv from ai
static void ocp_qp_daqp_unpack_vec(int n, struct blasfeo_dvec *sv, int ai, double *v)
{
    for (int i = 0; i < n; i++)
        v[i] = BLASFEO_DVECEL(sv, ai+i);
}



// Bit of the mask of ocp_qp_daqp_update_data: the weights of the soft constraints have changed
#define OCP_QP_DAQP_SOFT_WEIGHTS 64
// ... which constraints are soft has changed
#define OCP_QP_DAQP_SOFT_PATTERN 128

// v = e, setting *changed if v changes
static void ocp_qp_daqp_set(double *v, double e, int *changed)
{
    if (*v != e)
    {
        *v = e;
        *changed = 1;
    }
}

// Copy qp_in into the format of DAQP. Returns DAQP_MULTISTAGE_UPDATE_MATRICES if
// the matrices have changed (DAQP detects changes of the vectors itself), and
// OCP_QP_DAQP_SOFT_* for the soft constraints
static int ocp_qp_daqp_update_data(const ocp_qp_in *in, ocp_qp_daqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;
    DAQPMultistage *ocp = &mem->ocp;
    int N = dims->N;
    int *nx = dims->nx, *nu = dims->nu, *nb = dims->nb, *ng = dims->ng, *ns = dims->ns;
    int mat = 0, weights = 0, pattern = 0, c = 0;

    // the matrices of qp_in, used by DAQP directly (only their first nu+nx
    // rows: b and [r; q] are taken from the vectors b and rqz)
    ocp->BAbt = in->BAbt;
    ocp->RSQrq = in->RSQrq;
    ocp->DCt = in->DCt;
    for (int t = 0; t <= N; t++)
    {
        int nut = t < N ? nu[t] : 0, nux = nut + nx[t];
        // dynamics: BAbt = [B'; A'; b']
        if (t < N)
        {
            mat |= ocp_qp_daqp_sync_mat(nux, nx[t+1], in->BAbt+t, mem->cBAbt+t);
            ocp_qp_daqp_unpack_vec(nx[t+1], in->b+t, 0, ocp->b[t]);
        }
        // cost: RSQrq = [R S'; S Q; r' q'] (lower triangle), rqz = [r; q; z]
        mat |= ocp_qp_daqp_sync_mat(nux, nux, in->RSQrq+t, mem->cRSQrq+t);
        ocp_qp_daqp_unpack_vec(nut, in->rqz+t, 0, ocp->r[t]);
        ocp_qp_daqp_unpack_vec(nx[t], in->rqz+t, nut, ocp->q[t]);
        // general constraints: DCt = [D'; C']
        mat |= ocp_qp_daqp_sync_mat(nux, ng[t], in->DCt+t, mem->cDCt+t);
        // bounds: d = [lb; lg; -ub; -ug], with sides that are masked out unbounded
        const double *d = in->d[t].pa, *dm = in->d_mask[t].pa;
        int nbg = nb[t] + ng[t];
        double *lb = ocp->lb[t], *ub = ocp->ub[t], *lg = ocp->lg[t], *ug = ocp->ug[t];
        memcpy(ocp->idxb[t], in->idxb[t], nb[t]*sizeof(int));
        for (int i = 0; i < nb[t]; i++)
        {
            lb[i] = dm[i] != 0.0 ? d[i] : -DAQP_INF;
            ub[i] = dm[nbg+i] != 0.0 ? -d[nbg+i] : DAQP_INF;
        }
        for (int i = 0; i < ng[t]; i++)
        {
            lg[i] = dm[nb[t]+i] != 0.0 ? d[nb[t]+i] : -DAQP_INF;
            ug[i] = dm[nbg+nb[t]+i] != 0.0 ? -d[nbg+nb[t]+i] : DAQP_INF;
        }
        // soft constraints: slack j with the penalty 0.5*Z*s^2 + z*s, s >= s_min
        for (int i = 0; i < nbg; i++, c++)
        {
            int j = in->idxs_rev[t][i];
            int soft = j >= 0 ? DAQP_SOFT : 0;
            if (mem->sense[c] != soft)
            {
                mem->sense[c] = soft;
                pattern = 1;
            }
            if (!soft)
                continue;
            double Zl = MAX(1e-8, BLASFEO_DVECEL(in->Z+t, j));
            double Zu = MAX(1e-8, BLASFEO_DVECEL(in->Z+t, ns[t]+j));
            double zl = BLASFEO_DVECEL(in->rqz+t, nut+nx[t]+j);
            double zu = BLASFEO_DVECEL(in->rqz+t, nut+nx[t]+ns[t]+j);
            double sl = MAX(d[2*nbg+j], -zl/Zl);
            double su = MAX(d[2*nbg+ns[t]+j], -zu/Zu);
            double *lo = i < nb[t] ? ocp->lb[t]+i : ocp->lg[t]+i-nb[t];
            double *hi = i < nb[t] ? ocp->ub[t]+i : ocp->ug[t]+i-nb[t];
            if (*lo > -DAQP_INF) *lo -= sl;
            if (*hi < DAQP_INF) *hi += su;
            mem->s_l[c] = sl;
            mem->s_u[c] = su;
            ocp_qp_daqp_set(mem->rho_l+c, 1/Zl, &weights);
            ocp_qp_daqp_set(mem->rho_u+c, 1/Zu, &weights);
            ocp_qp_daqp_set(mem->w_l+c, MAX(0, zl + Zl*sl), &weights);
            ocp_qp_daqp_set(mem->w_u+c, MAX(0, zu + Zu*su), &weights);
        }
    }
    int soft = (weights ? OCP_QP_DAQP_SOFT_WEIGHTS : 0) | (pattern ? OCP_QP_DAQP_SOFT_PATTERN : 0);

    return (mat ? DAQP_MULTISTAGE_UPDATE_MATRICES : 0) | soft;
}



static void ocp_qp_daqp_fill_output(const ocp_qp_in *in, ocp_qp_out *out, ocp_qp_daqp_memory *mem)
{
    ocp_qp_dims *dims = in->dim;
    DAQPMultistage *ocp = &mem->ocp;
    int N = dims->N;
    int *nx = dims->nx, *nu = dims->nu, *nb = dims->nb, *ng = dims->ng, *ns = dims->ns;
    int zo = 0, co = 0, po = 0;

    for (int t = 0; t <= N; t++)
    {
        int nut = t < N ? nu[t] : 0, nbg = nb[t] + ng[t];
        // ux = [u; x; s_l; s_u] (z = [x; u] per stage)
        const double *x = mem->z+zo, *u = x+nx[t];
        blasfeo_pack_dvec(nut, mem->z+zo+nx[t], 1, out->ux+t, 0);
        blasfeo_pack_dvec(nx[t], mem->z+zo, 1, out->ux+t, nut);
        zo += nx[t] + nut;
        // lam = [lam_lb; lam_lg; lam_ub; lam_ug] (nonnegative). A side that is
        // masked out (unbounded) has no multiplier: a multiplier of rounding
        // size with that sign would otherwise multiply its huge slack in the
        // complementarity of acados
        blasfeo_dvecse(2*nbg, 0.0, out->lam+t, 0);
        for (int i = 0; i < nbg; i++)
        {
            double l = mem->lam[co+i];
            if (l < 0)
                BLASFEO_DVECEL(out->lam+t, i) = -l*BLASFEO_DVECEL(in->d_mask+t, i);
            else
                BLASFEO_DVECEL(out->lam+t, nbg+i) = l*BLASFEO_DVECEL(in->d_mask+t, nbg+i);
        }
        // slacks s = s_hat + t, and the multipliers of s >= s_min (stationarity
        // in s: Z*s + z - lam - lam_s = 0). The slack t of DAQP is nonzero only
        // beyond the linear weight w = z + Z*s_hat, where lam = w + Z*t, which
        // is exact in stationarity and complementarity (a violation within
        // the tolerance of DAQP is left in the constraint). A constraint that
        // DAQP cannot change (as one on a fixed x_0) has t = its violation
        if (ns[t] > 0)
            blasfeo_dvecse(2*ns[t], 0.0, out->lam+t, 2*nbg);
        for (int i = 0; i < nbg; i++)
        {
            int j = in->idxs_rev[t][i], c = co+i;
            if (j < 0)
                continue;
            double sl = mem->s_l[c], su = mem->s_u[c], l = mem->lam[c];
            double v = 0, lo, hi;
            if (i < nb[t])
            {
                int k = ocp->idxb[t][i];
                v = k < nut ? u[k] : x[k-nut];
                lo = ocp->lb[t][i]; hi = ocp->ub[t][i];
            }
            else
            {
                int k = i-nb[t];  // (DCt = [D'; C'])
                for (int r = 0; r < nut; r++) v += BLASFEO_DMATEL(in->DCt+t, r, k)*u[r];
                for (int r = 0; r < nx[t]; r++) v += BLASFEO_DMATEL(in->DCt+t, nut+r, k)*x[r];
                lo = ocp->lg[t][k]; hi = ocp->ug[t][k];
            }
            if (mem->setup_done && daqp_multistage_fixed_row(mem->work, c))
            {
                sl += MAX(0, lo - v);
                su += MAX(0, v - hi);
            }
            else if (l < 0 && -l > mem->w_l[c])
                sl += MAX(0, lo - v);
            else if (l > 0 && l > mem->w_u[c])
                su += MAX(0, v - hi);
            double Zl = 1/mem->rho_l[c], Zu = 1/mem->rho_u[c];
            double zl = BLASFEO_DVECEL(in->rqz+t, nut+nx[t]+j), zu = BLASFEO_DVECEL(in->rqz+t, nut+nx[t]+ns[t]+j);
            BLASFEO_DVECEL(out->ux+t, nut+nx[t]+j) = sl;
            BLASFEO_DVECEL(out->ux+t, nut+nx[t]+ns[t]+j) = su;
            BLASFEO_DVECEL(out->lam+t, 2*nbg+j) = Zl*sl + zl - BLASFEO_DVECEL(out->lam+t, i);
            BLASFEO_DVECEL(out->lam+t, 2*nbg+ns[t]+j) = Zu*su + zu - BLASFEO_DVECEL(out->lam+t, nbg+i);
        }
        co += nbg;
        if (t < N)
        {
            blasfeo_pack_dvec(nx[t+1], mem->pi+po, 1, out->pi+t, 0);
            po += nx[t+1];
        }
    }
}



void ocp_qp_daqp_terminate(void *config_, void *mem_, void *work_)
{
    ocp_qp_daqp_free_work(mem_);
}



int ocp_qp_daqp(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, void *work_)
{
    ocp_qp_in *qp_in = qp_in_;
    ocp_qp_out *qp_out = qp_out_;
    ocp_qp_daqp_opts *opts = opts_;
    ocp_qp_daqp_memory *mem = mem_;
    ocp_qp_dims *dims = qp_in->dim;
    DAQPWorkspace *work = mem->work;
    DAQPResult res;
    int flag;

    qp_info *info = (qp_info *) qp_out->misc;
    acados_timer tot_timer, qp_timer, interface_timer;
    acados_tic(&tot_timer);

    // not supported: inputs at the last stage, and slacks that soften several constraints
    int unsupported = dims->nu[dims->N] > 0;
    for (int t = 0; t <= dims->N; t++)
    {
        int nbg = dims->nb[t] + dims->ng[t];
        for (int i = 0; i < nbg; i++)
            for (int k = i+1; k < nbg; k++)
                if (qp_in->idxs_rev[t][i] >= 0 && qp_in->idxs_rev[t][i] == qp_in->idxs_rev[t][k])
                    unsupported = 1;
    }
    if (unsupported)
    {
        printf("\nerror: ocp_qp_daqp: nu[N] > 0, and slacks of several constraints, are not supported\n");
        mem->status = ACADOS_UNKNOWN;
        return ACADOS_UNKNOWN;
    }

    // data
    acados_tic(&interface_timer);
    int mask = ocp_qp_daqp_update_data(qp_in, mem);
    info->interface_time = acados_toc(&interface_timer);

    acados_tic(&qp_timer);
    if (mask & OCP_QP_DAQP_SOFT_PATTERN)
        ocp_qp_daqp_free_work(mem);  // (which constraints are soft is fixed at setup)
    // factorization (the working set is kept, and is a warm start)
    work->settings = opts->daqp_opts;
    flag = 1;
    if (mem->setup_done)
    {
        if (mask & OCP_QP_DAQP_SOFT_WEIGHTS)
        {
            // (the state of the soft constraints in the working set belongs to
            // the previous weights: it is dropped by a refactorization)
            daqp_set_soft_weights(work, mem->rho_l, mem->rho_u, mem->w_l, mem->w_u);
            mask |= DAQP_MULTISTAGE_UPDATE_MATRICES;
        }
        // cold start: an empty working set (but with the equality constraints
        // that DAQP keeps in it). With a refactorization, the update activates
        // the working set that is marked; otherwise it is activated here
        if (!opts->warm_start)
            for (int i = 0; i < work->m; i++)
                if (!(work->sense[i] & DAQP_IMMUTABLE))
                    work->sense[i] &= ~(DAQP_ACTIVE | DAQP_LOWER | DAQP_SLACK_FIXED);
        flag = daqp_multistage_update(work, mask & ~(OCP_QP_DAQP_SOFT_WEIGHTS | OCP_QP_DAQP_SOFT_PATTERN));
        if (flag == DAQP_EXIT_UNSUPPORTED)
            ocp_qp_daqp_free_work(mem);  // x_0 is no longer fixed (or now is): a new setup
        else if (flag < 0)
            ocp_qp_daqp_free_work(mem);  // (the factorization failed: set up anew at the next call)
        else if (!opts->warm_start && !(mask & DAQP_MULTISTAGE_UPDATE_MATRICES))
        {
            reset_daqp_workspace(work);
            flag = daqp_activate_constraints(work);
        }
    }
    if (!mem->setup_done)
    {
        flag = setup_daqp_multistage_mem(&mem->ocp, work, mem->daqp_mem, mem->daqp_mem_size);
        mem->setup_done = flag >= 0;
        int has_soft = 0;  // (the weights are only set for soft constraints)
        for (int i = 0; i < work->m; i++)
            if (mem->sense[i] & DAQP_SOFT)
                has_soft = 1;
        if (mem->setup_done && has_soft && !daqp_set_soft_weights(work, mem->rho_l, mem->rho_u, mem->w_l, mem->w_u))
        {
            printf("\nerror: ocp_qp_daqp: DAQP is built without individual soft weights (SOFT_WEIGHTS)\n");
            flag = DAQP_EXIT_UNSUPPORTED;
        }
    }

    // solve
    if (flag >= 0)
    {
        res.x = mem->z;
        res.lam = mem->lam;
        daqp_multistage_solve(&res, work);
        flag = res.exitflag;
        mem->iter = res.iter;
    }
    else
        mem->iter = 0;
    mem->time_qp_solver_call = acados_toc(&qp_timer);

    // fill qp_out (the slacks t are computed when needed)
    if (flag > 0)
        daqp_multistage_duals(work, mem->lam, mem->pi); // complete multipliers and costates
    ocp_qp_daqp_fill_output(qp_in, qp_out, mem);
    info->t_computed = 0;

    info->solve_QP_time = mem->time_qp_solver_call;
    info->total_time = acados_toc(&tot_timer);
    info->num_iter = mem->iter;

    int acados_status;
    if (flag == DAQP_EXIT_OPTIMAL || flag == DAQP_EXIT_SOFT_OPTIMAL || flag == DAQP_EXIT_OPTIMAL_INEXACT)
        acados_status = ACADOS_SUCCESS;
    else if (flag == DAQP_EXIT_ITERLIMIT)
        acados_status = ACADOS_MAXITER;
    else if (flag == DAQP_EXIT_INFEASIBLE)
        acados_status = ACADOS_INFEASIBLE;
    else if (flag == DAQP_EXIT_UNBOUNDED)
        acados_status = ACADOS_UNBOUNDED;
    else
        acados_status = ACADOS_UNKNOWN;
    if (opts->print_level > 0)
        printf("ocp_qp_daqp: exitflag %d, %d iterations (update mask %d)\n", flag, mem->iter, mask);
    mem->status = acados_status;
    return acados_status;
}



void ocp_qp_daqp_eval_forw_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_daqp_eval_forw_sens: not implemented yet\n");
    exit(1);
}



void ocp_qp_daqp_eval_adj_sens(void *config_, void *qp_in, void *seed, void *qp_out, void *opts_, void *mem_, void *work_)
{
    printf("\nerror: ocp_qp_daqp_eval_adj_sens: not implemented yet\n");
    exit(1);
}



void ocp_qp_daqp_solver_get(void *config_, void *qp_in_, void *qp_out_, void *opts_, void *mem_, const char *field, int stage, void* value, int size1, int size2)
{
    printf("\nerror: ocp_qp_daqp_solver_get: not implemented yet\n");
    exit(1);
}



void ocp_qp_daqp_config_initialize_default(void *config_)
{
    qp_solver_config *config = config_;

    config->opts_calculate_size = &ocp_qp_daqp_opts_calculate_size;
    config->opts_assign = &ocp_qp_daqp_opts_assign;
    config->opts_initialize_default = &ocp_qp_daqp_opts_initialize_default;
    config->opts_update = &ocp_qp_daqp_opts_update;
    config->opts_set = &ocp_qp_daqp_opts_set;
    config->opts_get = &ocp_qp_daqp_opts_get;
    config->memory_calculate_size = &ocp_qp_daqp_memory_calculate_size;
    config->memory_assign = &ocp_qp_daqp_memory_assign;
    config->memory_get = &ocp_qp_daqp_memory_get;
    config->workspace_calculate_size = &ocp_qp_daqp_workspace_calculate_size;
    config->evaluate = &ocp_qp_daqp;
    config->terminate = &ocp_qp_daqp_terminate;
    config->eval_forw_sens = &ocp_qp_daqp_eval_forw_sens;
    config->eval_adj_sens = &ocp_qp_daqp_eval_adj_sens;
    config->memory_reset = &ocp_qp_daqp_memory_reset;
    config->solver_get = &ocp_qp_daqp_solver_get;
}
