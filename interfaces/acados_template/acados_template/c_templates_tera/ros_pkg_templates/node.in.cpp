#include "{{ ros_opts.package_name }}/node.h"

namespace {{ ros_opts.package_name }}
{

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
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
    this->setup_parameter_handlers();
    this->declare_parameters();
    this->load_parameters();
    this->log_parameters();
    param_callback_handle_ = this->add_on_set_parameters_callback(
        std::bind(&{{ ClassName }}::on_parameter_update, this, std::placeholders::_1)
    );

    // --- Subscriber ---
    {%- if ns %}
    {%- set state_topic = "/" ~ ros_opts.namespace ~ "/state" %}
    {%- else %}
    {%- set state_topic = "/state" %}
    {%- endif %}
    state_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::State>(
        "{{ state_topic }}", 10,
        std::bind(&{{ ClassName }}::state_callback, this, std::placeholders::_1));
    {%- if ns %}
    {%- set references_topic = "/" ~ ros_opts.namespace ~ "/references" %}
    {%- else %}
    {%- set references_topic = "/references" %}
    {%- endif %}
    references_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::References>(
        "{{ references_topic }}", 10,
        std::bind(&{{ ClassName }}::references_callback, this, std::placeholders::_1));
    {%- if ns %}
    {%- set parameters_topic = "/" ~ ros_opts.namespace ~ "/parameters" %}
    {%- else %}
    {%- set parameters_topic = "/parameters" %}
    {%- endif %}
    parameters_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::Parameters>(
        "{{ parameters_topic }}", 10,
        std::bind(&{{ ClassName }}::parameters_callback, this, std::placeholders::_1));


    // --- Publisher ---
    {%- if ns %}
    {%- set input_topic = "/" ~ ros_opts.namespace ~ "/control_input" %}
    {%- else %}
    {%- set input_topic = "/control_input" %}
    {%- endif %}
    control_input_pub_ = this->create_publisher<{{ ros_opts.package_name }}_interface::msg::ControlInput>(
        "{{ input_topic }}", 10);

    // --- Init solver ---
    this->initialize_solver();
    this->start_control_timer({{ solver_options.Tsim }});
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

    this->set_cost_weights();
    this->set_constraints();
    {%- if has_slack %}
    void set_slack_weights();
    {%- endif %}

    RCLCPP_INFO(this->get_logger(), "Acados solver initialized successfully.");
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
    if (status == ACADOS_SUCCESS) {
        this->get_input(u0_.data(), 0);
        this->publish_input(u0_, status);
        {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
        first_solve_ = false;
        this->prepare_rti_solve();
        {%- endif %}
    } else {
        this->publish_input(u0_default_, status);
        if (status == ACADOS_NAN_DETECTED) agilex_acados_reset(ocp_capsule_, 1);
        {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
        first_solve_ = true;
        {%- endif %}
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
    // Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
    parameter_handlers_["{{ ros_opts.package_name }}.constraints.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.constraints.{{ field }}, res);
        };
    {%- endif %}
    {%- endfor %}

    // Weights
    {%- for field, param in cost %}
    {%- if param and (field is starting_with('W')) %}
    parameter_handlers_["{{ ros_opts.package_name }}.weights.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.weights.{{ field }}, res);
        };
    {%- endif %}
    {%- endfor %}
    {%- if has_slack %}

    // Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) %}
    parameter_handlers_["{{ ros_opts.package_name }}.slacks.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.slacks.{{ field }}, res);
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Solver Options
    parameter_handlers_["{{ ros_opts.package_name }}.solver_options.Tsim"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->config_.solver_options.Tsim = p.as_double();
            try {
                this->start_control_timer(this->config_.solver_options.Tsim);
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
    this->declare_parameter("{{ ros_opts.package_name }}.weights.{{ field }}", std::vector<double>{
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
    this->declare_parameter("{{ ros_opts.package_name }}.slacks.{{ field }}", std::vector<double>{ {{- param | join(sep=', ') -}} });
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Solver Options
    this->declare_parameter("{{ ros_opts.package_name }}.solver_options.Tsim", {{ solver_options.Tsim }});
}

void {{ ClassName }}::load_parameters() {
    // Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
    get_and_check_array_param(this, "{{ ros_opts.package_name }}.constraints.{{ field }}", config_.constraints.{{ field }});
    {%- endif %}
    {%- endfor %}

    // Weights
    {%- for field, param in cost %}
    {%- if param and (field is starting_with('W')) %}
    get_and_check_array_param(this, "{{ ros_opts.package_name }}.weights.{{ field }}", config_.weights.{{ field }});
    {%- endif %}
    {%- endfor %}
    {%- if has_slack %}

    // Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) %}
    get_and_check_array_param(this, "{{ ros_opts.package_name }}.slacks.{{ field }}", config_.slacks.{{ field }});
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Solver Options
    this->get_parameter("{{ ros_opts.package_name }}.solver_options.Tsim", config_.solver_options.Tsim);
}

void {{ ClassName }}::log_parameters() {
    const int label_width = 25;
    std::stringstream ss;

    ss << "\n----- {{ model.name | upper }} MPC Configuration -----";
    // Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ field }}" << " = " << config_.constraints.{{ field }};
    {%- endif %}
    {%- endfor %}

    // Weights
    {%- for field, param in cost %}
    {%- if param and (field is starting_with('W')) %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ field }}" << " = " << config_.weights.{{ field }};
    {%- endif %}
    {%- endfor %}
    {%- if has_slack %}

    // Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ field }}" << " = " << config_.slacks.{{ field }};
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    
    // Solver Options
    ss << "\n" << std::left << std::setw(label_width) << "Tsim" << " = " << config_.solver_options.Tsim;
    ss << "\n" << "--------------------------------------";

    RCLCPP_DEBUG(this->get_logger(), "%s", ss.str().c_str());
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
        } else {
            result.reason = "Update for unknown parameter '%s' received.", param_name.c_str();
            result.successful = false;
        }
    }

    if (result.successful){
        this->set_constraints();
        this->set_cost_weights();
        {% if has_slack %}
        this->set_slack_weights();
        {% endif %}
        this->log_parameters();
    }
    return result;
}


// --- Helpers ---
void {{ ClassName }}::start_control_timer(double period_seconds) {
    if (period_seconds <= 0.0) period_seconds = 0.02;
    auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::duration<double>(period_seconds));
    control_timer_ = this->create_wall_timer(
        period,
        std::bind(&{{ ClassName }}::control_loop, this));
}


// --- Acados Helpers ---
{%- if solver_options.nlp_solver_type == 'SQP_RTI' %}
void {{ ClassName }}::warmstart_solver_states(double *x0) {
    for (int i = 1; i <= AGILEX_N; ++i) {
        ocp_nlp_out_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, ocp_nlp_in_, i, "x", x0);
    }
}

int {{ ClassName }}::prepare_rti_solve() {
    int phase = PREPARATION;
    ocp_nlp_sqp_rti_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "rti_phase", &phase);
    int status = {{ model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS && status != ACADOS_READY) {
        first_solve_ = true;
        RCLCPP_ERROR(this->get_logger(), "Solver failed at preperation phase: %d", status);
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

void {{ ClassName }}::set_cost_weights() {
    {%- if dims.ny_0 > 0 %}
    auto W_0 = diag_from_vec(config_.weights.W_0);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "W", W_0.data());
    {%- endif %}
    {%- if dims.ny > 0 %}

    auto W = diag_from_vec(config_.weights.W);
    for (int i = 1; i < {{ model.name | upper }}_N; i++) {
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "W", W.data());
    }
    {%- endif %}
    {%- if dims.ny_e > 0 %}

    auto W_e = diag_from_vec(config_.weights.W_e);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, {{ model.name | upper }}_N, "W", W_e.data());
    {%- endif %}
}

{%- if has_slack %}
void {{ ClassName }}::set_slack_weights() {
    {%- if dims.ns_0 > 0 %}
    // Initial Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_0')) %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "{{ field }}", config_.slacks.{{ field }}.data());
    {%- endif %}
    {%- endfor %}

    {%- endif %}
    {%- if dims.ns > 0 %}
    // Stage Slacks
    for (int i = 1; i < {{ model.name | upper }}_N; i++) {
        {%- for field, param in cost %}
        {%- set field_l = field | lower %}
        {%- if param and (field_l is starting_with('z')) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i,  "{{ field }}", config_.slacks.{{ field }}.data());
        {%- endif %}
        {%- endfor %}
    }
    {%- endif %}
    {%- if dims.ns_e > 0 %}

    // Terminal Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_e')) %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "{{ field }}", config_.slacks.{{ field }}.data());
    {%- endif %}
    {%- endfor %}
    {%- endif %}
}
{%- endif %}

void {{ ClassName }}::set_constraints() {
    {%- if dims.nh_0 > 0 or dims.nphi_0 > 0 or dims.nsh_0 > 0 or dims.nsphi_0 > 0 %}
    // Initial Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_0')) and ('bx_0' not in field) %}
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "{{ field }}", config_.constraints.{{ field }}.data());
    {%- endif %}
    {%- endfor %}

    {%- endif %}
    {%- if dims.nbu > 0 or dims.nbx > 0 or dims.ng > 0 or dims.nh > 0 or dims.nphi > 0 or dims.nsbx > 0 or dims.nsg > 0 or dims.nsh > 0 or dims.nsphi > 0 %}
    // Stage Constraints
    for (int i=1; i < {{ model.name | upper }}_N; i++) {
        {%- for field, param in constraints %}
        {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
        ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, i, "{{ field }}", config_.constraints.{{ field }}.data());
        {%- endif %}
        {%- endfor %}
    }
    {%- endif %}
    {%- if dims.nbx_e > 0 or dims.ng_e > 0 or dims.nh_e > 0 or dims.nphi_e > 0 or dims.nsbx_e > 0 or dims.nsg_e > 0 or dims.nsh_e > 0 or dims.nsphi_e > 0 %}

    // Terminal Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and (field is ending_with('_e')) %}
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, {{ model.name | upper }}_N, "{{ field }}", config_.constraints.{{ field }}.data());
    {%- endif %}
    {%- endfor %}
    {%- endif %}
}

} // namespace {{ ros_opts.package_name }}


// --- Main Funktion ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<{{ ros_opts.package_name }}::{{ ClassName }}>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}