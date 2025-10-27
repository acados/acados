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
{%- set control_sequence_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic ~ "_sequence" %}
{%- else %}
{%- set control_topic = "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.state_topic %}
{%- set references_topic = "/" ~ ros_opts.reference_topic %}
{%- set parameters_topic = "/" ~ ros_opts.parameters_topic %}
{%- set control_sequence_topic = "/" ~ ros_opts.control_topic ~ "_sequence" %}
{%- endif %}
{%- set has_slack = dims.ns > 0 or dims.ns_0 > 0 or dims.ns_e > 0 %}
{%- set use_multithreading = ros_opts.threads is defined and ros_opts.threads > 1 %}
{{ ClassName }}::{{ ClassName }}()
    : Node("{{ ros_opts.node_name }}"), control_timer_(nullptr)
{
    RCLCPP_INFO(this->get_logger(), "Initializing {{ ros_opts.node_name | replace(from="_", to=" ") | title }}...");

    // --- default values ---
    config_ = {{ ClassName }}Config();
    {%- if constraints.has_x0 %}
    current_x_ = { {{- constraints.lbx_0 | join(sep=', ') -}} };
    {%- else %}
    current_x_.fill(0.0);
    {%- endif %}
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

    {%- if use_multithreading %}
    // --- Multithreading ---
    timer_group_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
    services_group_ = this->create_callback_group(rclcpp::CallbackGroupType::MutuallyExclusive);
    auto cb_options = rclcpp::SubscriptionOptions();
    cb_options.callback_group = services_group_;
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
        std::bind(&{{ ClassName }}::state_callback, this, std::placeholders::_1)
        {%- if use_multithreading -%}, cb_options {%- endif -%});
    references_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::References>(
        "{{ references_topic }}", 10,
        std::bind(&{{ ClassName }}::references_callback, this, std::placeholders::_1)
        {%- if use_multithreading -%}, cb_options {%- endif -%});
    {%- if dims.np > 0 %}
    parameters_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::Parameters>(
        "{{ parameters_topic }}", 10,
        std::bind(&{{ ClassName }}::parameters_callback, this, std::placeholders::_1)
        {%- if use_multithreading -%}, cb_options {%- endif -%});
    {%- endif %}

    // --- Publisher ---
    control_pub_ = this->create_publisher<{{ ros_opts.package_name }}_interface::msg::Control>(
        "{{ control_topic }}", 10);
    {%- if ros_opts.publish_control_sequence %}
    control_sequence_pub_ = this->create_publisher<{{ ros_opts.package_name }}_interface::msg::ControlSequence>(
        "{{ control_sequence_topic }}", 10);
    {%- endif %}

    // --- Init solver ---
    this->initialize_solver();
    this->apply_all_parameters_to_solver();
    this->start_control_timer(config_.ts);
}

{{ ClassName }}::~{{ ClassName }}() {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
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
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    ocp_capsule_ = {{ model.name }}_acados_create_capsule();
    int status = {{ model.name }}_acados_create(ocp_capsule_);
    if (status) {
        RCLCPP_FATAL(this->get_logger(), "{{ model.name }}acados_create() failed with status %d.", status);
        rclcpp::shutdown();
    }

    ocp_nlp_in_ = {{ model.name }}_acados_get_nlp_in(ocp_capsule_);
    ocp_nlp_out_ = {{ model.name }}_acados_get_nlp_out(ocp_capsule_);
    ocp_nlp_sens_ = {{ model.name }}_acados_get_sens_out(ocp_capsule_);
    ocp_nlp_config_ = {{ model.name }}_acados_get_nlp_config(ocp_capsule_);
    ocp_nlp_opts_ = {{ model.name }}_acados_get_nlp_opts(ocp_capsule_);
    ocp_nlp_dims_ = {{ model.name }}_acados_get_nlp_dims(ocp_capsule_);

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
        {%- if use_multithreading %}
        std::scoped_lock lock(data_mutex_);
        {%- endif %}
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

    {
        {%- if use_multithreading %}
        std::lock_guard<std::recursive_mutex> lock(solver_mutex_);

        {%- endif %}
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
}

void {{ ClassName }}::solver_status_behaviour(int status) {
    // publish u0 also if the solver failed
    this->get_control(u0_.data(), 0);
    this->publish_control(u0_, status);

    {%- if ros_opts.publish_control_sequence %}
    // publish full control sequence
    for (size_t i = 0; i < {{ solver_options.N_horizon }}; ++i) {
        this->get_control(u_seq_[i].data(), static_cast<int>(i));
    }
    this->publish_control_sequence(u_seq_, status);
    {%- endif %}

    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    // prepare for next iteration
    if (status == ACADOS_SUCCESS) {
        first_solve_ = false;
        this->prepare_rti_solve();
    }
    else {
        first_solve_ = true;
        if (config_.verbose) {
            {{ model.name }}_acados_print_stats(ocp_capsule_);
        }
        RCLCPP_INFO(this->get_logger(), "Resetting acados solver memory...");
        {{ model.name }}_acados_reset(ocp_capsule_, 1);
    }
    {%- endif %}
}


// --- ROS Callbacks ---
void {{ ClassName }}::state_callback(const {{ ros_opts.package_name }}_interface::msg::State::SharedPtr msg) {
    if (std::any_of(msg->x.begin(), msg->x.end(), [](double val){ return !std::isfinite(val); })) {
        RCLCPP_WARN_THROTTLE(
            this->get_logger(), *this->get_clock(), 2000,
            "State callback received NaN/Inf in 'x'. Ignoring message.");
        return;
    }

    {%- if use_multithreading %}
    std::scoped_lock lock(data_mutex_);
    {%- endif %}
    current_x_ = msg->x;
    RCLCPP_DEBUG_STREAM_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "State callback: x=" << current_x_);
}

void {{ ClassName }}::references_callback(const {{ ros_opts.package_name }}_interface::msg::References::SharedPtr msg) {
    {%- if dims.ny_0 > 0 %}
    if (std::any_of(msg->yref_0.begin(), msg->yref_0.end(), [](double val){ return !std::isfinite(val); })) {
        RCLCPP_WARN_THROTTLE(
            this->get_logger(), *this->get_clock(), 2000,
            "Reference callback received NaN/Inf in 'yref_0'. Ignoring message.");
        return;
    }
    {%- endif %}
    {%- if dims.ny > 0 %}
    if (std::any_of(msg->yref.begin(), msg->yref.end(), [](double val){ return !std::isfinite(val); })) {
        RCLCPP_WARN_THROTTLE(
            this->get_logger(), *this->get_clock(), 2000,
            "Reference callback received NaN/Inf in 'yref'. Ignoring message.");
        return;
    }
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    if (std::any_of(msg->yref_e.begin(), msg->yref_e.end(), [](double val){ return !std::isfinite(val); })) {
        RCLCPP_WARN_THROTTLE(
            this->get_logger(), *this->get_clock(), 2000,
            "Reference callback received NaN/Inf in 'yref_e'. Ignoring message.");
        return;
    }
    {%- endif %}

    {%- if use_multithreading %}
    std::scoped_lock lock(data_mutex_);
    {%- endif %}
    {%- if dims.ny_0 > 0 %} current_yref_0_ = msg->yref_0; {% endif %}
    {%- if dims.ny > 0 %}  current_yref_  = msg->yref;    {% endif %}
    {%- if dims.ny_e > 0 %} current_yref_e_ = msg->yref_e; {% endif %}
    {%- if dims.ny_0 > 0 %} RCLCPP_DEBUG_STREAM_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Refs callback: yref_0=" << current_yref_0_); {% endif %}
    {%- if dims.ny > 0 %}  RCLCPP_DEBUG_STREAM_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Refs callback: yref="   << current_yref_);   {% endif %}
    {%- if dims.ny_e > 0 %} RCLCPP_DEBUG_STREAM_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Refs callback: yref_e=" << current_yref_e_);  {% endif %}
}
{%- if dims.np > 0 %}

void {{ ClassName }}::parameters_callback(const {{ ros_opts.package_name }}_interface::msg::Parameters::SharedPtr msg) {
    if (std::any_of(msg->p.begin(), msg->p.end(), [](double val){ return !std::isfinite(val); })) {
        RCLCPP_WARN_THROTTLE(
            this->get_logger(), *this->get_clock(), 2000,
            "Parameter callback received NaN/Inf in 'p'. Ignoring message.");
        return;
    }

    {%- if use_multithreading %}
    std::scoped_lock lock(data_mutex_);
    {%- endif %}
    current_p_ = msg->p;
    RCLCPP_DEBUG_STREAM_THROTTLE(this->get_logger(), *this->get_clock(), 1000, "Params callback: p=" << current_p_);
}
{%- endif %}


// --- ROS Publisher ---
void {{ ClassName }}::publish_control(
        const std::array<double, {{ model.name | upper }}_NU>& u0, 
        int status
) {
    auto control = std::make_unique<{{ ros_opts.package_name }}_interface::msg::Control>();
    control->header.stamp = this->get_clock()->now();
    control->status = status;
    control->u = u0;
    control_pub_->publish(std::move(control));
}
{%- if ros_opts.publish_control_sequence %}

void {{ ClassName }}::publish_control_sequence(
        const std::array<std::array<double, {{ model.name | upper }}_NU>, {{ solver_options.N_horizon }}>& u_sequence, 
        int status
) {
    auto control_sequence = std::make_unique<{{ ros_opts.package_name }}_interface::msg::ControlSequence>();
    control_sequence->header.stamp = this->get_clock()->now();
    control_sequence->status = status;
    
    for (size_t i = 0; i < {{ solver_options.N_horizon }}; ++i) {
        control_sequence->control_sequence[i].u = u_sequence[i];
    }
    control_sequence_pub_->publish(std::move(control_sequence));
}
{%- endif %}


// --- ROS Parameter ---
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
    parameter_handlers_["constraints.{{ field }}"] =
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
    parameter_handlers_["constraints.{{ field }}"] =
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
    parameter_handlers_["constraints.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_constraint<{{ constraint_size }}>(p, res, "{{ field }}", std::vector<int>{ {{- model.name | upper -}}_N});
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Weights
    {%- if dims.ny_0 > 0 %}
    parameter_handlers_["cost.W_0"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NY0>(p, res, "W", std::vector<int>{0});
        };
    {%- endif %}
    {%- if dims.ny > 0 %}
    parameter_handlers_["cost.W"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            auto stages = range(1, {{ model.name | upper }}_N);
            this->update_cost<{{ model.name | upper }}_NY>(p, res, "W", stages);
        };
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    parameter_handlers_["cost.W_e"] =
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
    parameter_handlers_["cost.{{ field }}"] =
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
    parameter_handlers_["cost.{{ field }}"] =
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
    parameter_handlers_["cost.{{ field }}"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->update_cost<{{ model.name | upper }}_NSN>(p, res, "{{ field }}", std::vector<int>{ {{- model.name | upper -}}_N});
        };
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- endif %}

    // Solver Options
    parameter_handlers_["solver_options.print_level"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult&) {
            {%- if use_multithreading %}
            std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
            {%- endif %}
            int print_level = p.as_int();
            ocp_nlp_solver_opts_set(ocp_nlp_config_, ocp_nlp_opts_, "print_level", &print_level);
        };

    // Ros Configs
    parameter_handlers_["ts"] =
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
    parameter_handlers_["verbose"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult&) {
            this->config_.verbose = p.as_bool();
        };
}

void {{ ClassName }}::declare_parameters() {
    // Constraints
    {%- for field, param in constraints %}
    {%- if param and ((field is starting_with('l')) or (field is starting_with('u'))) and ('bx_0' not in field) %}
    this->declare_parameter("constraints.{{ field }}", std::vector<double>{ {{- param | join(sep=', ') -}} });
    {%- endif %}
    {%- endfor %}

    // Weights
    {%- for field, param in cost %}
    {%- if param and (field is starting_with('W')) %}
    this->declare_parameter("cost.{{ field }}", std::vector<double>{
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
    this->declare_parameter("cost.{{ field }}", std::vector<double>{ {{- param | join(sep=', ') -}} });
    {%- endif %}
    {%- endfor %}
    {%- endif %}

    // Solver Options
    this->declare_parameter("solver_options.print_level", {{ 0 }});

    // Ros Configs
    this->declare_parameter("ts", {{ solver_options.Tsim }});
    this->declare_parameter("verbose", false);
}

void {{ ClassName }}::load_parameters() {
    this->get_parameter("ts", config_.ts);
    this->get_parameter("verbose", config_.verbose);
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
            result.reason = "Update for unknown parameter '" + param_name + "' received.";
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

    {
        {%- if use_multithreading %}
        std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
        {%- endif %}
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

    {
        {%- if use_multithreading %}
        std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
        {%- endif %}
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
        std::bind(&{{ ClassName }}::control_loop, this)
        {%- if use_multithreading %}, timer_group_{%- endif -%});
}


// --- Acados Solver ---
{%- if solver_options.nlp_solver_type == 'SQP_RTI' %}
int {{ ClassName }}::prepare_rti_solve() {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
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
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
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
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    int status = {{ model.name }}_acados_solve(ocp_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Solver failed with status: %d", status);
    }
    return status;
}


// --- Acados Getter ---
void {{ ClassName }}::get_control(double* u, int stage) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    ocp_nlp_out_get(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, stage, "u", u);
}

void {{ ClassName }}::get_state(double* x, int stage) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    ocp_nlp_out_get(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, stage, "x", x);
}


// --- Acados Setter ---
void {{ ClassName }}::set_x0(double* x0) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    int lbx_status = ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "lbx", x0);
    int ubx_status = ocp_nlp_constraints_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, ocp_nlp_out_, 0, "ubx", x0);
    check_acados_status("set_x0 (lbx)", 0, lbx_status);
    check_acados_status("set_x0 (ubx)", 0, ubx_status);
}

void {{ ClassName }}::set_yref(double* yref, int stage) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    int status = ocp_nlp_cost_model_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_in_, stage, "yref", yref);
    check_acados_status("set_yref (yref)", stage, status);
}

{%- if dims.ny_0 > 0 %}
void {{ ClassName }}::set_yref0(double* yref0) {
    this->set_yref(yref0, 0);
}
{%- endif %}
{%- if dims.ny > 0 %}

void {{ ClassName }}::set_yrefs(double* yref) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    for (int i = 1; i < {{ model.name | upper }}_N; i++) {
        this->set_yref(yref, i);
    }
}
{%- endif %}
{%- if dims.ny_e > 0 %}

void {{ ClassName }}::set_yref_e(double* yrefN) {
    this->set_yref(yrefN, {{ model.name | upper }}_N);
}
{%- endif %}
{%- if dims.np > 0 %}

void {{ ClassName }}::set_ocp_parameter(double* p, size_t np, int stage) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    int status = {{ model.name }}_acados_update_params(ocp_capsule_, stage, p, np);
    check_acados_status("set_ocp_parameter (p)", stage, status);
}

void {{ ClassName }}::set_ocp_parameters(double* p, size_t np) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    for (int i = 0; i <= {{ model.name | upper }}_N; i++) {
        this->set_ocp_parameter(p, np, i);
    }
}
{%- endif %}
{%- if solver_options.nlp_solver_type == 'SQP_RTI' %}

void {{ ClassName }}::warmstart_solver_states(double *x0) {
    {%- if use_multithreading %}
    std::lock_guard<std::recursive_mutex> lock(solver_mutex_);
    {%- endif %}
    for (int i = 1; i <= {{ model.name | upper }}_N; ++i) {
        ocp_nlp_out_set(ocp_nlp_config_, ocp_nlp_dims_, ocp_nlp_out_, ocp_nlp_in_, i, "x", x0);
    }
}
{%- endif %}


// --- Acados Helper ---
bool {{ ClassName }}::check_acados_status(const char* field, int stage, int status) {
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "%s failed at stage %d: %d", field, stage, status);
        return false;
    }
    return true;
}

} // namespace {{ ros_opts.package_name }}


// --- Main ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);

    // Suppress ROS timer logging
    rcutils_logging_set_logger_level("rcl", RCUTILS_LOG_SEVERITY_WARN);
    rcutils_logging_set_logger_level("rclcpp", RCUTILS_LOG_SEVERITY_WARN);

    auto node = std::make_shared<{{ ros_opts.package_name }}::{{ ClassName }}>();
{%- if use_multithreading %}
    RCLCPP_INFO(node->get_logger(), "Using MultiThreadedExecutor with %d threads.", {{ ros_opts.threads }});
    rclcpp::executors::MultiThreadedExecutor executor(rclcpp::ExecutorOptions(), {{ ros_opts.threads }});
    executor.add_node(node);
    executor.spin();
{%- else %}
    RCLCPP_INFO(node->get_logger(), "Using SingleThreadedExecutor.");
    rclcpp::spin(node);
{%- endif %}
    rclcpp::shutdown();
    return 0;
}
