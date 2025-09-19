#include "{{ ros_opts.package_name }}/node.h"

namespace {{ ros_opts.package_name }}
{

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
{%- set control_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.parameters_topic %}
{%- else %}
{%- set control_topic = "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.parameters_topic %}
{%- endif %}
{%- set has_slack = dims.ns > 0 or dims.ns_0 > 0 or dims.ns_e > 0 %}
{{ ClassName }}::{{ ClassName }}()
    : Node("{{ ros_opts.node_name }}")
{
    RCLCPP_INFO(this->get_logger(), "Initializing {{ ros_opts.node_name | replace(from="_", to=" ") | title }}...");

    // --- default values ---
    config_ = {{ ClassName }}Config();
    {%- if dims.ny_0 > 0 %}
    current_yref_0_ = { {{- cost.yref_0 | join(sep=', ') -}} };
    {%- endif %}
    {%- if dims.ny > 0 %}
    current_yref_ = { {{- cost.yref | join(sep=', ') -}} };
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    current_yref_e_ = { {{- cost.yref_e | join(sep=', ') -}} };
    {%- endif %}
    {%- if dims.np > 0 %}
    current_p_ = { {{- parameter_values | join(sep=', ') -}} };
    {%- endif %}

    // --- Parameters ---
    this->declare_parameters();
    this->setup_parameter_handlers();
    param_callback_handle_ = this->add_on_set_parameters_callback(
        std::bind(&{{ ClassName }}::on_parameter_update, this, std::placeholders::_1));
    this->load_parameters();

    // --- Subscriber ---
    state_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::State>(
        "{{ state_topic }}", 10,
        std::bind(&{{ ClassName }}::state_callback, this, std::placeholders::_1));
    references_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::References>(
        "{{ references_topic }}", 10,
        std::bind(&{{ ClassName }}::references_callback, this, std::placeholders::_1));
    {%- if dims.np > 0 %}
    parameters_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::Parameters>(
        "{{ parameters_topic }}", 10,
        std::bind(&{{ ClassName }}::parameters_callback, this, std::placeholders::_1));
    {%- endif %}

    // --- Publisher ---
    control_input_pub_ = this->create_publisher<{{ ros_opts.package_name }}_interface::msg::ControlInput>(
        "{{ control_topic }}", 10);

    // --- Init solver ---
    this->initialize_solver();
    this->apply_all_parameters_to_solver();
    this->start_control_timer(config_.ts);
}

{{ ClassName }}::~{{ ClassName }}() {
    RCLCPP_INFO(this->get_logger(), "Shutting down and freeing Acados solver memory.");
    if (ocp_capsule_) {
        int status = {{ model.name }}_acados_free(ocp_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ model.name }}_acados_free() returned status %d.", status);
        }
        status = {{ model.name }}_acados_free_capsule(ocp_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ model.name }}_acados_free_capsule() returned status %d.", status);
        }
    }
}


// --- Core Methods ---
void {{ ClassName }}::initialize_solver() {
    ocp_capsule_ = {{ model.name }}_acados_create_capsule();
    int status = {{ model.name }}_acados_create(ocp_capsule_);
    if (status) {
        RCLCPP_FATAL(this->get_logger(), "{{ model.name }}acados_create() failed with status %d.", status);
        rclcpp::shutdown();
    }

    ocp_nlp_config_ = {{ model.name }}_acados_get_nlp_config(ocp_capsule_);
    ocp_nlp_dims_ = {{ model.name }}_acados_get_nlp_dims(ocp_capsule_);
    ocp_nlp_in_ = {{ model.name }}_acados_get_nlp_in(ocp_capsule_);
    ocp_nlp_out_ = {{ model.name }}_acados_get_nlp_out(ocp_capsule_);
    ocp_nlp_opts_ = {{ model.name }}_acados_get_nlp_opts(ocp_capsule_);

    RCLCPP_INFO(this->get_logger(), "acados solver initialized successfully.");
}

void {{ ClassName }}::control_loop() {
    // TODO: check for received msgs first
    std::array<double, {{ model.name | upper }}_NX> x0{};
    {%- if dims.ny_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NY0> yref0{};
    {%- endif %}
    {%- if dims.ny > 0 %}
    std::array<double, {{ model.name | upper }}_NY> yref{};
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    std::array<double, {{ model.name | upper }}_NYN> yrefN{};
    {%- endif %}
    {%- if dims.np > 0 %}
    std::array<double, {{ model.name | upper }}_NP> p{};
    {%- endif %}

    {
        std::scoped_lock lock(data_mutex_);
        x0 = current_x_;
        {%- if dims.ny_0 > 0 %}
        yref0 = current_yref_0_;
        {%- endif %}
        {%- if dims.ny > 0 %}
        yref = current_yref_;
        {%- endif %}
        {%- if dims.ny_e > 0 %}
        yrefN = current_yref_e_;
        {%- endif %}
        {%- if dims.np > 0 %}
        p = current_p_;
        {%- endif %}
    }

    // Update solver
    this->set_x0(x0.data());

    {%- if dims.ny_0 > 0 %}
    this->set_yref0(yref0.data());
    {%- endif %}
    {%- if dims.ny > 0 %}
    this->set_yrefs(yref.data());
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    this->set_yref_e(yrefN.data());
    {%- endif %}
    {%- if dims.np > 0 %}
    this->set_ocp_parameters(p.data(), p.size());
    {%- endif %}

    // Solve OCP
    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    int status;
    if (first_solve_) {
        this->warmstart_solver_states(x0.data());
        status = this->ocp_solve();
    }
    else {
        status = this->feedback_rti_solve();
    }

    {%- else %}
    int status = this->ocp_solve();
    {%- endif %}
    this->solver_status_behaviour(status);
}

void {{ ClassName }}::solver_status_behaviour(int status) {
    // publish u0 also if the solver failed
    this->get_input(u0_.data(), 0);
    this->publish_input(u0_, status);

    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    // prepare for next iteration
    if (status == ACADOS_SUCCESS) {
        first_solve_ = false;
        this->prepare_rti_solve();
    }
    else {
        first_solve_ = true;
    }
    {%- endif %}

    // reset solver if nan is detected
    if (status == ACADOS_NAN_DETECTED) {
        {{ model.name }}_acados_reset(ocp_capsule_, 1);
    }
}


// --- ROS Callbacks ---
void {{ ClassName }}::state_callback(const {{ ros_opts.package_name }}_interface::msg::State::SharedPtr msg) {
    std::scoped_lock lock(data_mutex_);
    std::copy_n(msg->x.begin(), {{ model.name | upper }}_NX, current_x_.begin());
}

void {{ ClassName }}::references_callback(const {{ ros_opts.package_name }}_interface::msg::References::SharedPtr msg) {
    std::scoped_lock lock(data_mutex_);
    {%- if dims.ny_0 > 0 %}
    std::copy_n(msg->yref_0.begin(), {{ model.name | upper }}_NY0, current_yref_0_.begin());
    {%- endif %}
    {%- if dims.ny > 0 %}
    std::copy_n(msg->yref.begin(), {{ model.name | upper }}_NY, current_yref_.begin());
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    std::copy_n(msg->yref_e.begin(), {{ model.name | upper }}_NYN, current_yref_e_.begin());
    {%- endif %}
}
{%- if dims.np > 0 %}

void {{ ClassName }}::parameters_callback(const {{ ros_opts.package_name }}_interface::msg::Parameters::SharedPtr msg) {
    std::scoped_lock lock(data_mutex_);
    std::copy_n(msg->p.begin(), {{ model.name | upper }}_NP, current_p_.begin());
}
{%- endif %}


// --- ROS Publisher ---
void {{ ClassName }}::publish_input(const std::array<double, {{ model.name | upper }}_NU>& u0, int status) {
    auto control_input = std::make_unique<{{ ros_opts.package_name }}_interface::msg::ControlInput>();
    control_input->header.stamp = this->get_clock()->now();
    control_input->header.frame_id = "{{ ros_opts.node_name }}_control_input";
    control_input->status = status;
    std::copy_n(u0.begin(), {{ model.name | upper }}_NU, control_input->u.begin());
    control_input_pub_->publish(std::move(control_input));
}


// --- Parameter Handling Methods ---
void {{ ClassName }}::setup_parameter_handlers() {
    {%- if dims.nh_0 > 0 or dims.nphi_0 > 0 or dims.nsh_0 > 0 or dims.nsphi_0 > 0 %}
    // Initial Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_0')) and ('bx_0' not in field) %}
    {%- set suffix = "" %}
    {%- if "h" in field and "s" not in field %}
        {%- set suffix = "_NH0" %}
    {%- elif "phi" in field and "s" not in field %}
        {%- set suffix = "_NPHI0" %}
    {%- elif "sh" in field %}
        {%- set suffix = "_NSH0" %}
    {%- elif "sphi" in field %}
        {%- set suffix = "_NSPHI0" %}
    {%- endif %}
    {%- set constraint_size = model.name ~ suffix | upper %}
    parameter_handlers_["{{ ros_opts.package_name }}.constraints.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_constraint<{{ constraint_size }}>(p, res, "{{ field }}", std::vector<int>{0});
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- if dims.nbu > 0 or dims.nbx > 0 or dims.ng > 0 or dims.nh > 0 or dims.nphi > 0 or dims.nsbx > 0 or dims.nsg > 0 or dims.nsh > 0 or dims.nsphi > 0 %}

    // Stage Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
    {%- set suffix = "" %}
    {%- if "bx" in field and "s" not in field %}
        {%- set suffix = "_NBX" %}
    {%- elif "bu" in field and "s" not in field %}
        {%- set suffix = "_NBU" %}
    {%- elif "h" in field and "s" not in field %}
        {%- set suffix = "_NH" %}
    {%- elif "phi" in field and "s" not in field %}
        {%- set suffix = "_NPHI" %}
    {%- elif "g" in field and "s" not in field %}
        {%- set suffix = "_NG" %}
    {%- elif "sbx" in field %}
        {%- set suffix = "_NSBX" %}
    {%- elif "sbu" in field %}
        {%- set suffix = "_NSBU" %}
    {%- elif "sh" in field %}
        {%- set suffix = "_NSH" %}
    {%- elif "sphi" in field %}
        {%- set suffix = "_NSPHI" %}
    {%- elif "sg" in field %}
        {%- set suffix = "_NSG" %}
    {%- endif %}
    {%- set constraint_size = model.name ~ suffix | upper %}
    parameter_handlers_["{{ ros_opts.package_name }}.constraints.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            auto stages = range(1, {{ model.name | upper }}_N);
            this->update_constraint<{{ constraint_size }}>(p, res, "{{ field }}", stages);
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- if dims.nbx_e > 0 or dims.ng_e > 0 or dims.nh_e > 0 or dims.nphi_e > 0 or dims.nsbx_e > 0 or dims.nsg_e > 0 or dims.nsh_e > 0 or dims.nsphi_e > 0 %}

    // Terminal Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_e')) %}
    {%- set suffix = "" %}
    {%- if "bx" in field and "s" not in field %}
        {%- set suffix = "_NBXN" %}
    {%- elif "bh" in field and "s" not in field %}
        {%- set suffix = "_NHN" %}
    {%- elif "phi" in field and "s" not in field %}
        {%- set suffix = "_NPHIN" %}
    {%- elif "g" in field and "s" not in field %}
        {%- set suffix = "_NGN" %}
    {%- elif "sbx" in field %}
        {%- set suffix = "_NSBXN" %}
    {%- elif "sh" in field %}
        {%- set suffix = "_NSHN" %}
    {%- elif "sphi" in field %}
        {%- set suffix = "_NSPHIN" %}
    {%- elif "sg" in field %}
        {%- set suffix = "_NSGN" %}
    {%- endif %}
    {%- set constraint_size = model.name ~ suffix | upper %}
    parameter_handlers_["{{ ros_opts.package_name }}.constraints.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_constraint<{{ constraint_size }}>(p, res, "{{ field }}", std::vector<int>{ {{- model.name | upper -}}_N});
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Weights
    {%- if dims.ny_0 > 0 %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.W_0"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NY0>(p, res, "W", std::vector<int>{0});
        };
    {%- endif %}
    {%- if dims.ny > 0 %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.W"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            auto stages = range(1, {{ model.name | upper }}_N);
            this->update_cost<{{ model.name | upper }}_NY>(p, res, "W", stages);
        };
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.W_e"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NYN>(p, res, "W", std::vector<int>{ {{- model.name | upper -}}_N});
        };
    {%- endif %}
    {%- if has_slack %}

    {%- if dims.ns_0 > 0 %}
    // Initial Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_0')) %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NS0>(p, res, "{{ field }}", std::vector<int>{0});
        };
    {%- endif %}
    {%- endfor %}

    {%- endif %}
    {%- if dims.ns > 0 %}
    // Stage Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            auto stages = range(1, {{ model.name | upper }}_N);
            this->update_cost<{{ model.name | upper }}_NS>(p, res, "{{ field }}", stages);
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- if dims.ns_e > 0 %}

    // Terminal Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_e')) %}
    parameter_handlers_["{{ ros_opts.package_name }}.cost.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NSN>(p, res, "{{ field }}", std::vector<int>{ {{- model.name | upper -}}_N});
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- endif %}

    // Solver Options
    parameter_handlers_["{{ ros_opts.package_name }}.ts"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->set_period(p.as_double());
            // Restart timer with the new period
            if (!this->is_running()) {
                // Assumes that the node is not yet ready and that the timer will be started later
                return;
            }
            try {
                this->start_control_timer(this->config_.ts);
            } catch (const std::exception& e) {
                res.reason = "Failed to start control timer, while setting parameter '" + p.get_name() + "': " + e.what();
                res.successful = false;
            }
        };
}

void {{ ClassName }}::declare_parameters() {
    // Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
    this->declare_parameter("{{ ros_opts.package_name }}.constraints.{{ field }}", std::vector<double>{ {{- param | join(sep=', ') -}} });
    {%- endif %}
    {%- endfor %}

    // Weights
    {%- for field, param in cost %}
    {%- if param and (field is starting_with('W')) %}
    this->declare_parameter("{{ ros_opts.package_name }}.cost.{{ field }}", std::vector<double>{
        {%- set n_diag = param | length -%}
        {%- for i in range(end=n_diag) -%}
            {{- param[i][i] -}}
            {%- if not loop.last %}, {% endif -%}
        {%- endfor -%}
    });
    {%- endif %}
    {%- endfor %}
    {%- if has_slack %}

    // Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) %}
    this->declare_parameter("{{ ros_opts.package_name }}.cost.{{ field }}", std::vector<double>{ {{- param | join(sep=', ') -}} });
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Solver Options
    this->declare_parameter("{{ ros_opts.package_name }}.ts", {{ solver_options.Tsim }});
}

void {{ ClassName }}::load_parameters() {
    this->get_parameter("{{ ros_opts.package_name }}.ts", config_.ts);
}

void {{ ClassName }}::apply_all_parameters_to_solver() {
    if (!ocp_capsule_) {
        RCLCPP_WARN(this->get_logger(), "apply_all_parameters_to_solver() called before solver init.");
        return;
    }
    rcl_interfaces::msg::SetParametersResult res;
    res.successful = true;

    for (auto & kv : parameter_handlers_) {
        const auto & name = kv.first;
        if (!this->has_parameter(name)) continue;
        auto param = this->get_parameter(name);
        kv.second(param, res);
        if (!res.successful) {
            RCLCPP_ERROR(this->get_logger(),
                "Failed to apply initial parameter '%s': %s",
                name.c_str(), res.reason.c_str());
            // reset flag for next parameter
            res.successful = true;
            res.reason.clear();
        }
    }
}

rcl_interfaces::msg::SetParametersResult {{ ClassName }}::on_parameter_update(
    const std::vector<rclcpp::Parameter>& params
) {
    rcl_interfaces::msg::SetParametersResult result;
    result.successful = true;

    for (const auto& param : params) {
        auto& param_name = param.get_name();

        if (parameter_handlers_.count(param_name)) {
            parameter_handlers_.at(param_name)(param, result);
            if (!result.successful) break;
        } else {
            result.reason = "Update for unknown parameter '%s' received.", param_name.c_str();
            result.successful = false;
        }
    }
    return result;
}

template <size_t N>
void {{ ClassName }}::update_param_array(
    const rclcpp::Parameter& param,
    std::array<double, N>& destination_array,
    rcl_interfaces::msg::SetParametersResult& result
) {
    auto values = param.as_double_array();

    if (values.size() != N) {
        result.successful = false;
        result.reason = "Parameter '" + param.get_name() + "' has size " +
                        std::to_string(values.size()) + ", but expected is " + std::to_string(N) + ".";
        return;
    }

    std::copy_n(values.begin(), N, destination_array.begin());
}

template<size_t N>
void {{ ClassName }}::update_constraint(
    const rclcpp::Parameter& param,
    rcl_interfaces::msg::SetParametersResult& result,
    const char* field,
    const std::vector<int>& stages
) {
    auto values = param.as_double_array();

    if (values.size() != N) {
        result.successful = false;
        result.reason = "Constraint '" + std::string(param.get_name()) + "' has size " +
                      std::to_string(values.size()) + ", but expected is " + std::to_string(N) + ".";
        return;
    }

    std::array<double, N> vec{};
    std::copy_n(values.begin(), N, vec.begin());

    for (int stage : stages) {
        int status = ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, stage, field, vec.data());

        if (status != ACADOS_SUCCESS) {
            result.successful = false;
            result.reason = "Acados solver failed to set cost field '" + std::string(field) +
                          "' for stage " + std::to_string(stage) + " (error code: " + std::to_string(status) + ")";
            return;
        }
    }
}

template<size_t N>
void {{ ClassName }}::update_cost(
    const rclcpp::Parameter& param,
    rcl_interfaces::msg::SetParametersResult& result,
    const char* field,
    const std::vector<int>& stages
) {
    const auto values = param.as_double_array();

    if (values.size() != N) {
        result.successful = false;
        result.reason = "Cost '" + std::string(param.get_name()) + "' has size " +
                      std::to_string(values.size()) + ", but expected is " + std::to_string(N) + ".";
        return;
    }

    std::array<double, N> vec;
    std::array<double, N * N> mat{};
    std::copy_n(values.begin(), N, vec.begin());
    double* data_ptr = nullptr;

    const bool is_weight = std::strcmp(field, "W") == 0;
    if (is_weight) {
        mat = diag_from_vec(vec);
        data_ptr = mat.data();
        RCLCPP_INFO_STREAM(this->get_logger(), "update cost field '" << field << "' mat(flat) = " << mat);
    } else {
        data_ptr = vec.data();
        RCLCPP_INFO_STREAM(this->get_logger(), "update cost field '" << field << "' values = " << vec);
    }

    for (int stage : stages) {
        int status = ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, stage, field, data_ptr);
        if (status != ACADOS_SUCCESS) {
            result.successful = false;
            result.reason = "Acados solver failed to set cost field '" + std::string(field) +
                          "' for stage " + std::to_string(stage) + " (error code: " + std::to_string(status) + ")";
            return;
        }
    }
}

// --- Helpers ---
void {{ ClassName }}::set_period(double period_seconds) {
    if (config_.ts == period_seconds) {
        // Nothing to do
        return;
    }
    RCLCPP_INFO_STREAM(
        this->get_logger(), "update period 'Ts' = " << period_seconds << "s");
    this->config_.ts = period_seconds;
    // Check period validity
    if (this->config_.ts <= 0.0) {
        this->config_.ts = {{ solver_options.time_steps[0] }};
        RCLCPP_WARN(this->get_logger(),
            "Control period must be positive, defaulting to {{ solver_options.time_steps[0] }}s, the first time step of the OCP.");
    }
}

void {{ ClassName }}::start_control_timer(double period_seconds) {
    this->set_period(period_seconds);
    // if timer already exists, restart with new period
    if (this->is_running()) {
        RCLCPP_WARN(this->get_logger(), "Control timer already running, restarting...");
        control_timer_->cancel();
    }
    RCLCPP_INFO_STREAM(this->get_logger(), "Starting control loop with period " << period_seconds << "s.");
    // create timer
    auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::duration<double>(period_seconds));
    control_timer_ = this->create_wall_timer(
        period,
        std::bind(&{{ ClassName }}::control_loop, this));
}


// --- Acados Helpers ---
{%- if solver_options.nlp_solver_type == 'SQP_RTI' %}
void {{ ClassName }}::warmstart_solver_states(double *x0) {
    for (int i = 1; i <= {{ model.name | upper }}_N; ++i) {
        ocp_nlp_out_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, ocp_nlp_in_, i, "x", x0);
    }
}

int {{ ClassName }}::prepare_rti_solve() {
    int phase = PREPARATION;
    ocp_nlp_sqp_rti_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "rti_phase", &phase);
    int status = {{ model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS && status != ACADOS_READY) {
        first_solve_ = true;
        RCLCPP_ERROR(this->get_logger(), "Solver failed at preparation phase: %d", status);
    }
    return status;
}

int {{ ClassName }}::feedback_rti_solve() {
    int phase = FEEDBACK;
    ocp_nlp_sqp_rti_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "rti_phase", &phase);
    int status = {{ model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed at feedback phase: %d", status);
    }
    return status;
}

{%- endif %}
int {{ ClassName }}::ocp_solve() {
    int status = {{ model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed with status: %d", status);
    }
    return status;
}

void {{ ClassName }}::get_input(double* u, int stage) {
    ocp_nlp_out_get(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, stage, "u", u);
}

void {{ ClassName }}::get_state(double* x, int stage) {
    ocp_nlp_out_get(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, stage, "x", x);
}

void {{ ClassName }}::set_x0(double* x0) {
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "lbx", x0);
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "ubx", x0);
}

{%- if dims.ny_0 > 0 %}
void {{ ClassName }}::set_yref0(double* yref0) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "yref", yref0);
}

{%- endif %}
{%- if dims.ny > 0 %}

void {{ ClassName }}::set_yref(double* yref, int stage) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, stage, "yref", yref);
}

void {{ ClassName }}::set_yrefs(double* yref) {
    for (int i = 1; i < {{ model.name | upper }}_N; i++) {
        this->set_yref(yref, i);
    }
}
{%- endif %}
{%- if dims.ny_e > 0 %}

void {{ ClassName }}::set_yref_e(double* yrefN) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, {{ model.name | upper }}_N, "yref", yrefN);
}
{%- endif %}
{%- if dims.np > 0 %}

void {{ ClassName }}::set_ocp_parameter(double* p, size_t np, int stage) {
    {{ model.name }}_acados_update_params(ocp_capsule_, stage, p, np);
}

void {{ ClassName }}::set_ocp_parameters(double* p, size_t np) {
    for (int i = 0; i <= {{ model.name | upper }}_N; i++) {
        this->set_ocp_parameter(p, np, i);
    }
}
{%- endif %}

} // namespace {{ ros_opts.package_name }}


// --- Main ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<{{ ros_opts.package_name }}::{{ ClassName }}>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
