#include "{{ package.name }}/{{ ros.node_name }}.h"

namespace {{ package.name }}
{
{% set ClassName = ros.node_name | replace('_', ' ') | title | replace(' ', '') %}
{% set has_slacks = acados.slacks.has_init or acados.slacks.has_term or acados.slacks.has_stage %}
{{ ClassName }}::{{ ClassName }}()
    : Node("{{ ros.node_name }}")
{
    RCLCPP_INFO(this->get_logger(), "Initializing {{ ros.node_name | replace('_', ' ') | title }}...");

    // --- default values ---
    config_ = {{ ClassName }}Config();
    {% if acados.solver.warmstart_first %}
    first_solve_ = true;
    {% endif %}
    u0_default_ = {};
    current_x_ = { {{ acados.x0.value | join(', ') }} };
    {% if acados.references.yref_0.value %}
    current_yref_0_ = { {{ acados.references.yref_0.value | join(', ') }} };
    {% endif %}
    {% if acados.references.yref.value %}
    current_yref_ = { {{ acados.references.yref.value | join(', ') }} };
    {% endif %}
    {% if acados.references.yref_e.value %}
    current_yref_e_ = { {{ acados.references.yref_e.value | join(', ') }} };
    {% endif %}
    {% if acados.parameter_values.value %}
    current_p_ = { {{ acados.parameter_values.value | join(', ') }} };
    {% endif %}

    // --- Parameters ---
    this->setup_parameter_handlers();
    this->declare_parameters();
    this->load_parameters();
    this->log_parameters();
    param_callback_handle_ = this->add_on_set_parameters_callback(
        std::bind(&{{ ClassName }}::on_parameter_update, this, std::placeholders::_1)
    );

    // --- Subscriber ---
    {% for sub in ros.subscribers %}
        {% if sub.msg_type is not none and sub.msg_type != 'None' %}
    {{ (sub.name | lower | replace(' ', '_')) }}_sub_ = this->create_subscription<{{ cpp_type(sub.msg_type) }}>(
        "{{ sub.topic }}", 10,
        std::bind(&{{ ClassName }}::{{ sub.callback | default((sub.name ~ '_callback')) }}, this, std::placeholders::_1));
        {% endif %}
    {% endfor %}

    // --- Publisher ---
    {% for pub in ros.publishers %}
        {% if pub.msg_type is not none and pub.msg_type != 'None' %}
    {{ (pub.name | lower | replace(' ', '_')) }}_pub_ = this->create_publisher<{{ cpp_type(pub.msg_type) }}>(
        "{{ pub.topic }}", {{ pub.queue_size }});
        {% endif %}
    {% endfor %}
    {% if package.with_markers == true %}
    marker_pub_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
        "visualization_marker_array", 10);
    {% endif %}

    // --- Init solver ---
    this->initialize_solver();
    this->start_control_timer({{ acados.solver.Tsim }});
}

{{ ClassName }}::~{{ ClassName }}() {
    RCLCPP_INFO(this->get_logger(), "Shutting down and freeing Acados solver memory.");
    if (ocp_capsule_) {
        int status = {{ acados.model.name }}_acados_free(ocp_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ acados.model.name }}_acados_free() returned status %d.", status);
        }
        status = {{ acados.model.name }}_acados_free_capsule(ocp_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ acados.model.name }}_acados_free_capsule() returned status %d.", status);
        }
    }
}


// --- Core Methods ---
void {{ ClassName }}::initialize_solver() {
    ocp_capsule_ = {{ acados.model.name }}_acados_create_capsule();
    int status = {{ acados.model.name }}_acados_create(ocp_capsule_);
    if (status) {
        RCLCPP_FATAL(this->get_logger(), "{{ acados.model.name }}acados_create() failed with status %d.", status);
        rclcpp::shutdown();
    }

    ocp_nlp_config_ = {{ acados.model.name }}_acados_get_nlp_config(ocp_capsule_);
    ocp_nlp_dims_ = {{ acados.model.name }}_acados_get_nlp_dims(ocp_capsule_);
    ocp_nlp_in_ = {{ acados.model.name }}_acados_get_nlp_in(ocp_capsule_);
    ocp_nlp_out_ = {{ acados.model.name }}_acados_get_nlp_out(ocp_capsule_);
    ocp_nlp_opts_ = {{ acados.model.name }}_acados_get_nlp_opts(ocp_capsule_);

    this->set_cost_weights();
    this->set_constraints();
    {% if has_slacks %}
    void set_slack_weights();
    {% endif %}

    RCLCPP_INFO(this->get_logger(), "Acados solver initialized successfully.");
}

void {{ ClassName }}::control_loop() {
    // TODO: check for received msgs first
    std::array<double, {{ acados.model.name | upper }}_NX> x0{}; 
    {% if acados.references.yref_0.value %}
    std::array<double, {{ acados.model.name | upper }}_NY0> yref0{}; 
    {% endif %}
    {% if acados.references.yref.value %}
    std::array<double, {{ acados.model.name | upper }}_NY> yref{};
    {% endif %}
    {% if acados.references.yref_e.value %}
    std::array<double, {{ acados.model.name | upper }}_NYN> yrefN{};
    {% endif %}
    {% if acados.parameter_values.value %}
    std::array<double, {{ acados.model.name | upper }}_NP> p{};
    {% endif %}

    {
        std::scoped_lock lock(data_mutex_);
        x0 = current_x_;
        {% if acados.references.yref_0.value %}
        yref0 = current_yref_0_;
        {% endif %}
        {% if acados.references.yref.value %}
        yref = current_yref_;
        {% endif %}
        {% if acados.references.yref_e.value %}
        yrefN = current_yref_e_;
        {% endif %}
        {% if acados.parameter_values.value %}
        p = current_p_;
        {% endif %}
    } 
    
    // Update solver
    this->set_x0(x0.data());
    {% if acados.references.yref_0.value %}
    this->set_yref0(yref0.data());
    {% endif %}
    {% if acados.references.yref.value %}
    this->set_yrefs(yref.data());
    {% endif %}
    {% if acados.references.yref_e.value %}
    this->set_yref_e(yrefN.data());
    {% endif %}
    {% if acados.parameter_values.value %}
    this->set_ocp_parameters(p.data(), p.size());
    {% endif %}

    {% if acados.solver.warmstart_first %}
    if (first_solve_) {
        this->warmstart_states(x0.data());
        first_solve_ = false;
    }
    {% endif %}
    {% if acados.solver.warmstart %}
    this->warmstart_states(x0.data());
    {% endif %}

    // Solve OCP
    {% if acados.solver.nlp_solver_type == "SQP_RTI" %} 
    int status = this->feedback_rti_solve();
    {% else %}
    int status = this->ocp_solve();
    {% endif %}
    if (status == ACADOS_SUCCESS) {
        std::array<double, {{ acados.model.name | upper }}_NU> u0;
        this->get_input(u0.data(), 0);
        this->publish_input(u0);
        {% if package.with_markers == true %}
        visualize_markers();
        {% endif %}
    } else {
        this->publish_input(u0_default_);
        RCLCPP_INFO(this->get_logger(), "Publishing default input.");
    }
    {% if acados.solver.nlp_solver_type == "SQP_RTI" %}

    this->prepare_rti_solve();
    {% endif %}
}


// --- ROS Callbacks ---
{% for sub in ros.subscribers %}
    {% if sub.msg_type is not none and sub.msg_type != 'None' %}
void {{ ClassName }}::{{ sub.callback | default((sub.name ~ '_callback')) }}(const {{ cpp_type(sub.msg_type) }}::SharedPtr msg) {
    std::scoped_lock lock(data_mutex_);
    // TODO: make a copy of all relevant data to call in the controll loop
}
    {% endif %}
{% endfor %}


// --- ROS Publisher ---
void {{ ClassName }}::publish_input(const std::array<double, {{ acados.model.name | upper }}_NU>& u0) {
    // TODO: publish the input with the correct message
    // auto cmd_vel = std::make_unique<geometry_msgs::msg::Twist>();
    // cmd_vel->linear.x = u0[0];
    // cmd_vel->angular.z = u0[1];
    // cmd_vel_pub_->publish(std::move(cmd_vel));
}


// --- Parameter Handling Methods ---
void {{ ClassName }}::setup_parameter_handlers() {
    // Constraints
    {% for field, param in acados.constraints.items() %}
    {% if param.value %}
    parameter_handlers_["{{ package.name }}.constraints.{{ param.name }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.constraints.{{ param.name }}, res);
        };
    {% endif %}
    {% endfor %}

    // Weights
    {% for field, param in acados.weights.items() %}
    {% if param.value %}
    parameter_handlers_["{{ package.name }}.weights.{{ param.name }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.weights.{{ param.name }}, res);
        };
    {% endif %}
    {% endfor %}
    {% if has_slacks %}

    // Slacks
    {% for field, param in acados.slacks.items() %}
    {% if param.value %}
    parameter_handlers_["{{ package.name }}.slacks.{{ param.name }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            update_param_array(p, this->config_.slacks.{{ param.name }}, res);
        };
    {% endif %}
    {% endfor %}
    {% endif %}

    // Solver Options
    parameter_handlers_["{{ package.name }}.solver.Tsim"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->config_.solver_options.Tsim = p.as_double();
            try {
                this->start_control_timer(this->config_.solver_options.Tsim);
            } catch (const std::exception& e) {
                res.reason = "Failed to start control timer, while setting parameter '" + p.get_name() + "': " + e.what();
                res.successful = false;
            }
        };

    // Other Parameters
    {% for param in ros.parameters %}
    parameter_handlers_["{{ param.name }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            try {
                {% if param.type == "string" %}
                this->config_.{{ param.name | snake }} = p.as_string();
                {% elif param.type == "double" or param.type == "float" %}
                this->config_.{{ param.name | snake }} = p.as_double();
                {% elif param.type == "int" %}
                this->config_.{{ param.name | snake }} = p.as_int();
                {% elif param.type == "bool" %}
                this->config_.{{ param.name | snake }} = p.as_bool();
                {% endif %}
            } catch (const std::exception& e) {
                res.reason = "Failed to set parameter '" + p.get_name() + "': " + e.what();
                res.successful = false;
            }
        };
    {% endfor %}
}

void {{ ClassName }}::declare_parameters() {
    // Constraints
    {% for field, param in acados.constraints.items() %}
        {% if param.value %}
    this->declare_parameter("{{ package.name }}.constraints.{{ param.name }}", std::vector<double>{ {{ param.value | join(', ') }} });
        {% endif %}
    {% endfor %}

    // Weights
    {% for field, param in acados.weights.items() %}
        {% if param.value %}
    this->declare_parameter("{{ package.name }}.weights.{{ param.name }}", std::vector<double>{ {{ param.value | join(', ') }} });
        {% endif %}
    {% endfor %}

    // Solver Options
    this->declare_parameter("{{ package.name }}.solver.Tsim", {{ acados.solver.Tsim }});

    // Other Parameters
    {% for param in ros.parameters %}
        {% if param.type == "string" %}
    this->declare_parameter("{{ param.name }}", "{{ param.value }}");
        {% else %}
    this->declare_parameter("{{ param.name }}", {{ param.value }});
        {% endif %}
    {% endfor %}
}

void {{ ClassName }}::load_parameters() {
    // Constraints
    {% for field, param in acados.constraints.items() %}
    {% if param.value %}
    get_and_check_array_param(this, "{{ package.name }}.constraints.{{ param.name }}", config_.constraints.{{ param.name }});
    {% endif %}
    {% endfor %}

    // Weights
    {% for field, param in acados.weights.items() %}
    {% if param.value %}
    get_and_check_array_param(this, "{{ package.name }}.weights.{{ param.name }}", config_.weights.{{ param.name }});
    {% endif %}
    {% endfor %}

    // Solver Options
    this->get_parameter("{{ package.name }}.solver.Tsim", config_.solver_options.Tsim);

    // Other Parameters
    {% for param in ros.parameters %}
    this->get_parameter("{{ package.name }}.{{ param.name }}", config_.{{ param.name | snake }});
    {% endfor %}
}

void {{ ClassName }}::log_parameters() {
    const int label_width = 25;
    std::stringstream ss;

    // Bauen Sie den gesamten String zusammen
    ss << "\n----- {{ acados.model.name | upper }} MPC Configuration -----";
    // Constraints
    {% for field, param in acados.constraints.items() %}
    {% if param.value %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ param.log_label }}" << " = " << config_.constraints.{{ param.name }};
    {% endif %}
    {% endfor %}

    // Weights
    {% for field, param in acados.weights.items() %}
    {% if param.value %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ param.log_label }}" << " = " << config_.weights.{{ param.name }};
    {% endif %}
    {% endfor %}
    
    // Solver Options
    ss << "\n" << std::left << std::setw(label_width) << "Tsim" << " = " << config_.solver_options.Tsim;

    // Other Parameters
    {% for param in ros.parameters %}
    ss << "\n" << std::left << std::setw(label_width) << "{{ param.name | replace('_', ' ') | title }}" << " = " << config_.{{ param.name | snake }};
    {% endfor %}
    ss << "\n" << "--------------------------------------";

    // Geben Sie den finalen String auf einmal aus
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
        {% if has_slacks %}
        void set_slack_weights();
        {% endif %}
        this->log_parameters();
    }
    return result;
}
{% if package.with_markers == true %}


// --- ROS Visualizer ---
void {{ ClassName }}::visualize_markers() {
    // TODO: collect the correct markers
    // auto marker_points = get_marker_points(
    //     some_geometry_points,
    //     "map",
    //     this->get_clock(),
    //     Color{0.0, 1.0, 0.0},
    //     "something",
    //     0,
    //     0.05,
    //     1.0,
    //     visualization_msgs::msg::Marker::LINE_STRIP
    // );

    // publish_marker_array(marker_pub_, marker_points);
}
{% endif %}


// --- Helpers ---
void {{ ClassName }}::start_control_timer(double rate_hz) {
    if (rate_hz <= 0.0) rate_hz = 50.0;
    auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::duration<double>(1.0 / rate_hz));
    control_timer_ = this->create_wall_timer(
        period,
        std::bind(&{{ ClassName }}::control_loop, this)
    );
}


// --- Acados Helpers ---
{% if acados.solver.nlp_solver_type == 'SQP_RTI' %}
int {{ ClassName }}::prepare_rti_solve() {
    int phase = PREPARATION;
    ocp_nlp_sqp_rti_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "rti_phase", &phase);
    int status = {{ acados.model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS && status != ACADOS_READY) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed at preperation phase: %d", status);
    }
    return status;
}

int {{ ClassName }}::feedback_rti_solve() {
    int phase = FEEDBACK;
    ocp_nlp_sqp_rti_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "rti_phase", &phase);
    int status = {{ acados.model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed at feedback phase: %d", status);
    }
    return status;
}
{% else %}
int {{ ClassName }}::ocp_solve() {
    int status = {{ acados.model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed with status: %d", status);
    }
    return status;
}
{% endif %}

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

void {{ ClassName }}::set_yref0(double* yref0) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "yref", yref0);
}

void {{ ClassName }}::set_yref(double* yref, int stage) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, stage, "yref", yref);
}

void {{ ClassName }}::set_yref_e(double* yrefN) {
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, {{ acados.model.name | upper }}_N, "yref", yrefN);
}

void {{ ClassName }}::set_yrefs(double* yref) {
    for (int i = 0; i <= {{ acados.model.name | upper }}_N; i++) {
        this->set_yref(yref, i);
    }
}
{% if acados.parameter_values.value %}

void {{ ClassName }}::set_ocp_parameter(double* p, size_t np, int stage) {
    {{ acados.model.name }}_acados_update_params(ocp_capsule_, stage, p, np);
}

void {{ ClassName }}::set_ocp_parameters(double* p, size_t np) {
    for (int i = 0; i <= {{ acados.model.name | upper }}_N; i++) {
        this->set_ocp_parameter(p, np, i);
    }
}
{% endif %}

void {{ ClassName }}::set_cost_weights() {
    {% if acados.weights.W_0.value %}
    auto W_0 = diag_from_vec(config_.weights.W_0);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "W", W_0.data());
    {% endif %}
    {% if acados.weights.W.value %}

    auto W = diag_from_vec(config_.weights.W);
    for (int i = 1; i < {{ acados.model.name | upper }}_N; i++) {
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "W", W.data());
    }
    {% endif %}
    {% if acados.weights.W_e.value %}

    auto W_e = diag_from_vec(config_.weights.W_e);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, {{ acados.model.name | upper }}_N, "W", W_e.data());
    {% endif %}
    return;
}

{% if has_slacks %}
void {{ ClassName }}::set_slack_weights() {
    {% if acados.slacks.Zl_0.value %}
    auto Zl_0 = diag_from_vec(config_.slacks.Zl_0);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "Zl", Zl_0.data());
    {% endif %}
    {% if acados.slacks.Zu_0.value %}
    auto Zu_0 = diag_from_vec(config_.slacks.Zu_0);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "Zu", Zu_0.data());
    {% endif %}
    {% if acados.slacks.zl_0.value %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "zl",  config_.slacks.zl_0.data());
    {% endif %}
    {% if acados.slacks.zu_0.value %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "zu", config_.slacks.zu_0.data());
    {% endif %}
    {% if acados.slacks.Zl.value or acados.slacks.Zu.value or acados.slacks.zl.value or acados.slacks.zu.value %}

    {% if acados.slacks.Zl.value %}
    auto Zl = diag_from_vec(config_.slacks.Zl);
    {% endif %}
    {% if acados.slacks.Zu.value %}
    auto Zu = diag_from_vec(config_.slacks.Zu);
    {% endif %}
    for (int i = 1; i < {{ acados.model.name | upper }}_N; i++) {
        {% if acados.slacks.Zl.value %}
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "Zl", Zl.data());
        {% endif %}
        {% if acados.slacks.Zu.value %}
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "Zu", Zu.data());
        {% endif %}
        {% if acados.slacks.zl.value %}
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "zl", config_.slacks.zl.data());
        {% endif %}
        {% if acados.slacks.zu.value %}
        ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, i, "zu", config_.slacks.zu.data());
        {% endif %}
    }
    {% endif %}
    {% if acados.slacks.Zl_e.value or acados.slacks.Zu_e.value or acados.slacks.zl_e.value or acados.slacks.zu_e.value %}

    {% if acados.slacks.Zl_e.value %}
    auto Zl_e = diag_from_vec(config_.slacks.Zl_e);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "Zl", Zl_e.data());
    {% endif %}
    {% if acados.slacks.Zu_e.value %}
    auto Zu_e = diag_from_vec(config_.slacks.Zu_e);
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "Zu", Zu_e.data());
    {% endif %}
    {% if acados.slacks.zl_e.value %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "zl",  config_.slacks.zl_e.data());
    {% endif %}
    {% if acados.slacks.zu_e.value %}
    ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, 0, "zu", config_.slacks.zu_e.data());
    {% endif %}
    {% endif %}
    return;
}
{% endif %}

void {{ ClassName }}::set_constraints() {
    {% if acados.constraints.has_init %}
    // Initial Constraints
    {% for field, param in acados.constraints.items() %}
    {% if param.value and param.name.endswith('_0') %}
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "{{ param.name }}", config_.constraints.{{ param.name }}.data());
    {% endif %}
    {% endfor %}
    {% endif %}
    {% if acados.constraints.has_stage %}

    // Stage Constraints
    for (int i=1; i < {{ acados.model.name | upper }}_N; i++) {
        {% for field, param in acados.constraints.items() %}
        {% if param.value and (not param.name.endswith('_0')) and (not param.name.endswith('_e')) %}
        ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, i, "{{ param.name }}", config_.constraints.{{ param.name }}.data());
        {% endif %}
        {% endfor %}
    }
    {% endif %}
    {% if acados.constraints.has_term %}

    // Terminal Constraints
    {% for field, param in acados.constraints.items() %}
    {% if param.value and param.name.endswith('_e') %}
    ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, {{ acados.model.name | upper }}_N, "{{ param.name }}", config_.constraints.{{ param.name }}.data());
    {% endif %}
    {% endfor %}
    {% endif %}
    return;
}
{% if acados.solver.warmstart or acados.solver.warmstart_first %}

void {{ ClassName }}::warmstart_states(double* x0) {
    for (int i = 1; i <= {{ acados.model.name | upper }}_N; ++i) {
        ocp_nlp_out_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, ocp_nlp_in_, i, "x", x0);
    }
}

void {{ ClassName }}::warmstart_inputs(double* u0) {
    for (int i = 0; i < {{ acados.model.name | upper }}_N; ++i) {
        ocp_nlp_out_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, ocp_nlp_in_, i, "u", u0);
    }
}
{% endif %}

} // namespace {{ package.name }}


// --- Main Funktion ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<{{ package.name }}::{{ ClassName }}>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}