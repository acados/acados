#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "{{ ros_opts.package_name }}/node.h"

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{%- set NodeType = ros_opts.package_name ~ "::" ~ ClassName %}
{%- set ConfigType = ros_opts.package_name ~ "::" ~ ClassName ~ "Config" %}
{%- set has_slack = dims.ns > 0 or dims.ns_0 > 0 or dims.ns_e > 0 %}
{%- set has_W_0 = dims.ny_0 > 0 and ((cost.cost_type_0 == "LINEAR_LS") or (cost.cost_type_0 == "NONLINEAR_LS")) %}
{%- set has_W = dims.ny > 0 and ((cost.cost_type == "LINEAR_LS") or (cost.cost_type == "NONLINEAR_LS")) %}
{%- set has_W_e = dims.ny_e > 0 and ((cost.cost_type_e == "LINEAR_LS") or (cost.cost_type_e == "NONLINEAR_LS")) %}

namespace
{

template <typename Tag>
struct MemberAccessor {
    static inline typename Tag::type ptr = nullptr;
};

template <typename Tag, typename Tag::type M>
struct AccessGrant {
    struct Setter {
        Setter() { MemberAccessor<Tag>::ptr = M; }
    };
    static Setter setter;
};

template <typename Tag, typename Tag::type M>
typename AccessGrant<Tag, M>::Setter AccessGrant<Tag, M>::setter;

struct NlpConfigTag { using type = ocp_nlp_config* {{ NodeType }}::*; };
struct NlpDimsTag   { using type = ocp_nlp_dims* {{ NodeType }}::*; };
struct NlpInTag     { using type = ocp_nlp_in* {{ NodeType }}::*; };
struct ConfigTag    { using type = {{ ConfigType }} {{ NodeType }}::*; };

template struct AccessGrant<NlpConfigTag, &{{ NodeType }}::ocp_nlp_config_>;
template struct AccessGrant<NlpDimsTag,   &{{ NodeType }}::ocp_nlp_dims_>;
template struct AccessGrant<NlpInTag,     &{{ NodeType }}::ocp_nlp_in_>;
template struct AccessGrant<ConfigTag,    &{{ NodeType }}::config_>;

[[maybe_unused]] std::vector<double> perturb_bound(const std::vector<double>& orig, bool is_lower) {
    std::vector<double> out = orig;
    const double sign = is_lower ? -1.0 : 1.0;
    for (size_t i = 0; i < out.size(); ++i) {
        const double delta = std::abs(out[i]) * 0.1 + 1.0 + 0.1 * static_cast<double>(i);
        out[i] += sign * delta;
    }
    return out;
}

[[maybe_unused]] void expect_constraint_matches(
    ocp_nlp_config* cfg,
    ocp_nlp_dims* dims,
    ocp_nlp_in* in,
    const std::vector<int>& stages,
    const char* acados_field,
    const std::vector<double>& expected
) {
    std::vector<double> actual(expected.size(), 0.0);
    for (int stage : stages) {
        ocp_nlp_constraints_model_get(cfg, dims, in, stage, acados_field, actual.data());
        for (size_t i = 0; i < expected.size(); ++i) {
            EXPECT_DOUBLE_EQ(actual[i], expected[i])
                << "Constraint mismatch for '" << acados_field << "' at stage " << stage << ", index " << i;
        }
    }
}

[[maybe_unused]] void expect_weight_col_major_matches(
    ocp_nlp_config* cfg,
    ocp_nlp_dims* dims,
    ocp_nlp_in* in,
    const std::vector<int>& stages,
    size_t dim,
    const std::vector<double>& expected_row_major
) {
    ASSERT_EQ(expected_row_major.size(), dim * dim);
    std::vector<double> actual_col_major(dim * dim, 0.0);
    for (int stage : stages) {
        ASSERT_EQ(ocp_nlp_cost_model_get(cfg, dims, in, stage, "W", actual_col_major.data()), ACADOS_SUCCESS);
        for (size_t r = 0; r < dim; ++r) {
            for (size_t c = 0; c < dim; ++c) {
                EXPECT_DOUBLE_EQ(actual_col_major[c * dim + r], expected_row_major[r * dim + c])
                    << "Weight W mismatch at stage " << stage << " (" << r << ", " << c << ")";
            }
        }
    }
}

[[maybe_unused]] void expect_cost_vec_matches(
    ocp_nlp_config* cfg,
    ocp_nlp_dims* dims,
    ocp_nlp_in* in,
    const std::vector<int>& stages,
    const char* acados_field,
    const std::vector<double>& expected
) {
    std::vector<double> actual(expected.size(), 0.0);
    for (int stage : stages) {
        ASSERT_EQ(ocp_nlp_cost_model_get(cfg, dims, in, stage, acados_field, actual.data()), ACADOS_SUCCESS);
        for (size_t i = 0; i < expected.size(); ++i) {
            EXPECT_DOUBLE_EQ(actual[i], expected[i])
                << "Cost mismatch for '" << acados_field << "' at stage " << stage << ", index " << i;
        }
    }
}

[[maybe_unused]] std::vector<int> make_stage_range(int start_stage, int end_exclusive) {
    std::vector<int> stages;
    for (int s = start_stage; s < end_exclusive; ++s) {
        stages.push_back(s);
    }
    return stages;
}

class OcpNodeSolverInspectionTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        if (!rclcpp::ok()) {
            rclcpp::init(0, nullptr);
        }
    }

    static void TearDownTestSuite() {
        if (rclcpp::ok()) {
            rclcpp::shutdown();
        }
    }

    void SetUp() override {
        // Construct without spinning an executor so control_timer_ never fires concurrently
        // and set_parameter() invokes on_parameter_update() synchronously in this thread.
        node_ = std::make_shared<{{ NodeType }}>();
    }

    void TearDown() override {
        node_.reset();
    }

    ocp_nlp_config* nlp_config() const { return node_.get()->*MemberAccessor<NlpConfigTag>::ptr; }
    ocp_nlp_dims* nlp_dims() const { return node_.get()->*MemberAccessor<NlpDimsTag>::ptr; }
    ocp_nlp_in* nlp_in() const { return node_.get()->*MemberAccessor<NlpInTag>::ptr; }
    const {{ ConfigType }}& node_config() const { return node_.get()->*MemberAccessor<ConfigTag>::ptr; }

    std::shared_ptr<{{ NodeType }}> node_;
};

TEST_F(OcpNodeSolverInspectionTest, SolverPointersAndHelperUtilsValid) {
    ASSERT_NE(nlp_config(), nullptr);
    ASSERT_NE(nlp_dims(), nullptr);
    ASSERT_NE(nlp_in(), nullptr);

    const std::array<double, 4> row_major{1.0, 2.0, 3.0, 4.0};
    const auto col_major = ::row_major_to_col_major<double, 2>(row_major);
    const std::array<double, 4> expected_col_major{1.0, 3.0, 2.0, 4.0};
    EXPECT_EQ(col_major, expected_col_major);
}

TEST_F(OcpNodeSolverInspectionTest, CostWeightsAndSlacksStoredAccuratelyAcrossStages) {
    {%- if has_W_0 %}
    {
        const size_t dim = {{ name | upper }}_NY0;
        const std::vector<int> stages{0};
        const std::vector<double> init_w{
        {%- for row in cost.W_0 -%}
            {%- for val in row -%}{{ val }}{% if not loop.last %}, {% endif %}{%- endfor -%}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
        };
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, init_w);

        std::vector<double> full_mod = init_w;
        for (size_t r = 0; r < dim; ++r) {
            for (size_t c = 0; c < dim; ++c) {
                full_mod[r * dim + c] += (r == c) ? (1.5 + 0.1 * static_cast<double>(r)) : 0.25;
            }
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W_0", full_mod)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, full_mod);

        std::vector<double> diag_vals(dim);
        std::vector<double> diag_expanded(dim * dim, 0.0);
        for (size_t i = 0; i < dim; ++i) {
            diag_vals[i] = init_w[i * dim + i] + 2.5 + static_cast<double>(i);
            diag_expanded[i * dim + i] = diag_vals[i];
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W_0", diag_vals)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);

        std::vector<double> bad_len(dim * dim + 1, 1.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.W_0", bad_len)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);
    }
    {%- endif %}

    {%- if has_W %}
    {
        const size_t dim = {{ name | upper }}_NY;
        const std::vector<int> stages = make_stage_range(1, {{ name | upper }}_N);
        const std::vector<double> init_w{
        {%- for row in cost.W -%}
            {%- for val in row -%}{{ val }}{% if not loop.last %}, {% endif %}{%- endfor -%}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
        };
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, init_w);

        std::vector<double> full_mod = init_w;
        for (size_t r = 0; r < dim; ++r) {
            for (size_t c = 0; c < dim; ++c) {
                full_mod[r * dim + c] += (r == c) ? (1.5 + 0.1 * static_cast<double>(r)) : 0.25;
            }
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W", full_mod)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, full_mod);

        std::vector<double> diag_vals(dim);
        std::vector<double> diag_expanded(dim * dim, 0.0);
        for (size_t i = 0; i < dim; ++i) {
            diag_vals[i] = init_w[i * dim + i] + 2.5 + static_cast<double>(i);
            diag_expanded[i * dim + i] = diag_vals[i];
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W", diag_vals)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);

        std::vector<double> bad_len(dim * dim + 1, 1.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.W", bad_len)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);
    }
    {%- endif %}

    {%- if has_W_e %}
    {
        const size_t dim = {{ name | upper }}_NYN;
        const std::vector<int> stages{ {{ name | upper }}_N };
        const std::vector<double> init_w{
        {%- for row in cost.W_e -%}
            {%- for val in row -%}{{ val }}{% if not loop.last %}, {% endif %}{%- endfor -%}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
        };
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, init_w);

        std::vector<double> full_mod = init_w;
        for (size_t r = 0; r < dim; ++r) {
            for (size_t c = 0; c < dim; ++c) {
                full_mod[r * dim + c] += (r == c) ? (1.5 + 0.1 * static_cast<double>(r)) : 0.25;
            }
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W_e", full_mod)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, full_mod);

        std::vector<double> diag_vals(dim);
        std::vector<double> diag_expanded(dim * dim, 0.0);
        for (size_t i = 0; i < dim; ++i) {
            diag_vals[i] = init_w[i * dim + i] + 2.5 + static_cast<double>(i);
            diag_expanded[i * dim + i] = diag_vals[i];
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.W_e", diag_vals)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);

        std::vector<double> bad_len(dim * dim + 1, 1.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.W_e", bad_len)).successful);
        expect_weight_col_major_matches(nlp_config(), nlp_dims(), nlp_in(), stages, dim, diag_expanded);
    }
    {%- endif %}

    {%- if has_slack %}
    {%- if dims.ns_0 > 0 %}
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_0')) %}
    {
        const std::vector<int> stages{0};
        const std::vector<double> init_z{ {{- param | join(sep=', ') -}} };
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_0', to='') }}", init_z);

        std::vector<double> mod_z = init_z;
        for (size_t i = 0; i < mod_z.size(); ++i) {
            mod_z[i] += 3.5 + 0.5 * static_cast<double>(i);
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", mod_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_0', to='') }}", mod_z);

        std::vector<double> bad_z(mod_z.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", bad_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_0', to='') }}", mod_z);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    {%- if dims.ns > 0 %}
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
    {
        const std::vector<int> stages = make_stage_range(1, {{ name | upper }}_N);
        const std::vector<double> init_z{ {{- param | join(sep=', ') -}} };
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", init_z);

        std::vector<double> mod_z = init_z;
        for (size_t i = 0; i < mod_z.size(); ++i) {
            mod_z[i] += 3.5 + 0.5 * static_cast<double>(i);
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", mod_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", mod_z);

        std::vector<double> bad_z(mod_z.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", bad_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", mod_z);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    {%- if dims.ns_e > 0 %}
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_e')) %}
    {
        const std::vector<int> stages{ {{ name | upper }}_N };
        const std::vector<double> init_z{ {{- param | join(sep=', ') -}} };
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_e', to='') }}", init_z);

        std::vector<double> mod_z = init_z;
        for (size_t i = 0; i < mod_z.size(); ++i) {
            mod_z[i] += 3.5 + 0.5 * static_cast<double>(i);
        }
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", mod_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_e', to='') }}", mod_z);

        std::vector<double> bad_z(mod_z.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("cost.{{ field }}", bad_z)).successful);
        expect_cost_vec_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field | replace(from='_e', to='') }}", mod_z);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- endif %}
}

TEST_F(OcpNodeSolverInspectionTest, ConstraintsUpdatedAndRejectedAtomicallyAcrossStages) {
    // Initial non-slack constraints (stage 0)
    {%- if (constraints.has_x0 == false and dims.nbx_0 > 0) or dims.nh_0 > 0 or dims.nphi_0 > 0 %}
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_0')) and (('bx_0' not in field) or (constraints.has_x0 == false)) and ('s' not in field) %}
    {
        const std::vector<int> stages{0};
        const std::vector<double> init_val{ {{- param | join(sep=', ') -}} };
        const char* acados_field = "{{ field | replace(from='_0', to='') }}";
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, init_val);

        const auto mod_val = perturb_bound(init_val, {% if field is starting_with('l') %}true{% else %}false{% endif %});
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", mod_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, mod_val);

        std::vector<double> bad_val(mod_val.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", bad_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, mod_val);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Stage non-slack constraints (bu & g cover stages 0..N-1; bx, h, phi cover stages 1..N-1)
    {%- if dims.nbu > 0 or dims.nbx > 0 or dims.ng > 0 or dims.nh > 0 or dims.nphi > 0 %}
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is not ending_with('_0')) and (field is not ending_with('_e')) and ('s' not in field) %}
    {%- set start_stage = 1 %}
    {%- if ("bu" in field) or ("g" in field) %}
        {%- set start_stage = 0 %}
    {%- endif %}
    {
        const std::vector<int> stages = make_stage_range({{ start_stage }}, {{ name | upper }}_N);
        const std::vector<double> init_val{ {{- param | join(sep=', ') -}} };
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", init_val);

        const auto mod_val = perturb_bound(init_val, {% if field is starting_with('l') %}true{% else %}false{% endif %});
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", mod_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", mod_val);

        std::vector<double> bad_val(mod_val.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", bad_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, "{{ field }}", mod_val);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Terminal non-slack constraints (stage N)
    {%- if dims.nbx_e > 0 or dims.ng_e > 0 or dims.nh_e > 0 or dims.nphi_e > 0 %}
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_e')) and ('s' not in field) %}
    {
        const std::vector<int> stages{ {{ name | upper }}_N };
        const std::vector<double> init_val{ {{- param | join(sep=', ') -}} };
        const char* acados_field = "{{ field | replace(from='_e', to='') }}";
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, init_val);

        const auto mod_val = perturb_bound(init_val, {% if field is starting_with('l') %}true{% else %}false{% endif %});
        EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", mod_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, mod_val);

        std::vector<double> bad_val(mod_val.size() + 1, 0.0);
        EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("constraints.{{ field }}", bad_val)).successful);
        expect_constraint_matches(nlp_config(), nlp_dims(), nlp_in(), stages, acados_field, mod_val);
    }
    {%- endif %}
    {%- endfor %}
    {%- endif %}
}

TEST_F(OcpNodeSolverInspectionTest, RosConfigAndOptionsEnforceValidation) {
    const double default_ts = static_cast<double>({{ solver_options.Tsim }});
    EXPECT_DOUBLE_EQ(node_config().ts, default_ts);
    EXPECT_EQ(node_config().threads, {{ ros_opts.threads | default(value=1) }});
    EXPECT_FALSE(node_config().verbose);

    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("ts", default_ts)).successful);
    EXPECT_DOUBLE_EQ(node_config().ts, default_ts);

    const double new_ts = default_ts * 2.0;
    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("ts", new_ts)).successful);
    EXPECT_DOUBLE_EQ(node_config().ts, new_ts);

    EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("ts", 0.0)).successful);
    EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("ts", -0.5)).successful);
    EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("ts", std::numeric_limits<double>::quiet_NaN())).successful);
    EXPECT_FALSE(node_->set_parameter(rclcpp::Parameter("ts", std::numeric_limits<double>::infinity())).successful);
    EXPECT_DOUBLE_EQ(node_config().ts, new_ts);

    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("verbose", true)).successful);
    EXPECT_TRUE(node_config().verbose);
    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("verbose", false)).successful);
    EXPECT_FALSE(node_config().verbose);

    const int default_print_level = {{ solver_options.print_level | default(value=0) }};
    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("solver_options.print_level", default_print_level == 0 ? 1 : 0)).successful);
    EXPECT_TRUE(node_->set_parameter(rclcpp::Parameter("solver_options.print_level", default_print_level)).successful);

    EXPECT_FALSE(node_->has_parameter("nonexistent_parameter"));
    EXPECT_THROW(
        node_->set_parameter(rclcpp::Parameter("nonexistent_parameter", 1.0)),
        rclcpp::exceptions::ParameterNotDeclaredException
    );
    EXPECT_THROW(
        node_->set_parameters({rclcpp::Parameter("nonexistent_parameter", 1.0)}),
        rclcpp::exceptions::ParameterNotDeclaredException
    );
}

} // namespace
