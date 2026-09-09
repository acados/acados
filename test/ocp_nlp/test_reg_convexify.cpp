/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */

#include <cmath>
#include <cstdlib>

#include "acados/ocp_nlp/ocp_nlp_reg_convexify.h"
#include "blasfeo_d_aux.h"
#include "blasfeo_d_aux_ext_dep.h"
#include "catch/include/catch.hpp"


namespace
{

class ConvexifyFixture
{
 public:
    ConvexifyFixture()
    {
        config_memory_ = std::malloc(ocp_nlp_reg_config_calculate_size());
        config_ = static_cast<ocp_nlp_reg_config *>(ocp_nlp_reg_config_assign(config_memory_));
        ocp_nlp_reg_convexify_config_initialize_default(config_);

        constexpr int N = 1;
        dims_memory_ = std::malloc(config_->dims_calculate_size(N));
        dims_ = config_->dims_assign(N, dims_memory_);
        dims_->nx[0] = 1;
        dims_->nx[1] = 1;
        dims_->nu[0] = 2;
        dims_->nu[1] = 0;

        opts_memory_ = std::malloc(config_->opts_calculate_size());
        opts_ = config_->opts_assign(opts_memory_);
        config_->opts_initialize_default(config_, dims_, opts_);

        memory_memory_ = std::malloc(config_->memory_calculate_size(config_, dims_, opts_));
        memory_ = config_->memory_assign(config_, dims_, opts_, memory_memory_);

        blasfeo_allocate_dmat(4, 3, &RSQrq_[0]);
        blasfeo_allocate_dmat(2, 1, &RSQrq_[1]);
        blasfeo_allocate_dmat(4, 1, &BAbt_[0]);
        blasfeo_allocate_dvec(3, &rq_[0]);
        blasfeo_allocate_dvec(1, &rq_[1]);
        blasfeo_allocate_dvec(1, &b_[0]);

        blasfeo_dgese(4, 3, 0.0, &RSQrq_[0], 0, 0);
        blasfeo_dgese(2, 1, 0.0, &RSQrq_[1], 0, 0);
        blasfeo_dgese(4, 1, 0.0, &BAbt_[0], 0, 0);
        blasfeo_dvecse(3, 0.0, &rq_[0], 0);
        blasfeo_dvecse(1, 0.0, &rq_[1], 0);
        blasfeo_dvecse(1, 0.0, &b_[0], 0);

        // Construct a control Hessian that becomes indefinite after the lower-triangular update.
        // Convexify should detect and regularize it. R = [[0.25, -0.75], [-0.75, 0.25]]
        const double delta = static_cast<ocp_nlp_reg_convexify_opts *>(opts_)->delta;
        BLASFEO_DMATEL(&RSQrq_[0], 0, 0) = 1.0;
        BLASFEO_DMATEL(&RSQrq_[0], 1, 1) = 1.0;
        BLASFEO_DMATEL(&RSQrq_[0], 2, 2) = 1.0;
        BLASFEO_DMATEL(&RSQrq_[1], 0, 0) = delta - 0.75;
        BLASFEO_DMATEL(&BAbt_[0], 0, 0) = 1.0;
        BLASFEO_DMATEL(&BAbt_[0], 1, 0) = 1.0;

        config_->memory_set_RSQrq_ptr(dims_, RSQrq_, memory_);
        config_->memory_set_rq_ptr(dims_, rq_, memory_);
        config_->memory_set_BAbt_ptr(dims_, BAbt_, memory_);
        config_->memory_set_b_ptr(dims_, b_, memory_);
    }

    ~ConvexifyFixture()
    {
        blasfeo_free_dvec(&b_[0]);
        blasfeo_free_dvec(&rq_[1]);
        blasfeo_free_dvec(&rq_[0]);
        blasfeo_free_dmat(&BAbt_[0]);
        blasfeo_free_dmat(&RSQrq_[1]);
        blasfeo_free_dmat(&RSQrq_[0]);
        std::free(memory_memory_);
        std::free(opts_memory_);
        std::free(dims_memory_);
        std::free(config_memory_);
    }

    void regularize() { config_->regularize(config_, dims_, opts_, memory_); }

    void regularize_lhs() { config_->regularize_lhs(config_, dims_, opts_, memory_); }

    void regularize_rhs() { config_->regularize_rhs(config_, dims_, opts_, memory_); }

    void check_control_hessian(const char *phase) const
    {
        const double r00 = BLASFEO_DMATEL(&RSQrq_[0], 0, 0);
        const double r01 = BLASFEO_DMATEL(&RSQrq_[0], 0, 1);
        const double r10 = BLASFEO_DMATEL(&RSQrq_[0], 1, 0);
        const double r11 = BLASFEO_DMATEL(&RSQrq_[0], 1, 1);

        INFO(phase);
        CHECK(std::isfinite(r00));
        CHECK(std::isfinite(r01));
        CHECK(std::isfinite(r10));
        CHECK(std::isfinite(r11));
        CHECK(r01 == Approx(r10).margin(1e-12));
        CHECK(r00 > 0.0);
        CHECK(r11 > 0.0);
        CHECK(r00 * r11 - r10 * r10 > 0.0);
    }

 private:
    void *config_memory_;
    void *dims_memory_;
    void *opts_memory_;
    void *memory_memory_;
    ocp_nlp_reg_config *config_;
    ocp_nlp_reg_dims *dims_;
    void *opts_;
    void *memory_;
    struct blasfeo_dmat RSQrq_[2];
    struct blasfeo_dmat BAbt_[1];
    struct blasfeo_dvec rq_[2];
    struct blasfeo_dvec b_[1];
};

}  // namespace


TEST_CASE("convexify uses the updated lower triangle", "[convexify]")
{
    ConvexifyFixture fixture;

    SECTION("full regularization")
    {
        fixture.regularize();
        fixture.check_control_hessian("after regularize");
    }

    SECTION("split regularization")
    {
        fixture.regularize_lhs();
        fixture.check_control_hessian("after regularize_lhs");
        fixture.regularize_rhs();
        fixture.check_control_hessian("after regularize_rhs");
    }
}
