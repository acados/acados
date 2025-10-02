#include "{{ ros_opts.package_name }}/node.h"

namespace {{ ros_opts.package_name }}
{

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{%- set ns = ros_opts.namespace | lower | trim(chars='/') | replace(from=" ", to="_") %}
{%- if ns %}
{%- set control_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.namespace ~ "/" ~ ros_opts.state_topic %}
{%- else %}
{%- set control_topic = "/" ~ ros_opts.control_topic %}
{%- set state_topic = "/" ~ ros_opts.state_topic %}
{%- endif %}
{{ ClassName }}::{{ ClassName }}()
    : Node("{{ ros_opts.node_name }}")
{
    RCLCPP_INFO(this->get_logger(), "Initializing {{ ros_opts.node_name | replace(from="_", to=" ") | title }}...");

    // --- default values ---
    config_ = {{ ClassName }}Config();

    // --- Parameters ---
    this->declare_parameters();
    this->setup_parameter_handlers();
    param_callback_handle_ = this->add_on_set_parameters_callback(
        std::bind(&{{ ClassName }}::on_parameter_update, this, std::placeholders::_1));
    this->load_parameters();

    // --- Subscriber ---
    control_sub_ = this->create_subscription<{{ ros_opts.package_name }}_interface::msg::ControlInput>(
        "{{ control_topic }}", 10,
        std::bind(&{{ ClassName }}::control_callback, this, std::placeholders::_1));

    // --- Publisher ---
    state_pub_ = this->create_publisher<{{ ros_opts.package_name }}_interface::msg::State>(
        "{{ state_topic }}", 10);

    // --- Init simulator ---
    this->initialize_simulator();
    this->start_integration_timer(config_.ts);
}

{{ ClassName }}::~{{ ClassName }}() {
    RCLCPP_INFO(this->get_logger(), "Shutting down and freeing Acados simulator memory.");
    if (sim_capsule_) {
        int status = {{ model.name }}_acados_sim_free(sim_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ model.name }}_acados_sim_free() returned status %d.", status);
        }
        status = {{ model.name }}_acados_sim_solver_free_capsule(sim_capsule_);
        if (status) {
            RCLCPP_ERROR(this->get_logger(), "{{ model.name }}_acados_sim_solver_free_capsule() returned status %d.", status);
        }
    }
}


// --- Core Methods ---
void {{ ClassName }}::initialize_simulator() {
    sim_capsule_ = {{ model.name }}_acados_sim_solver_create_capsule();
    int status = {{ model.name }}_acados_sim_create(sim_capsule_);
    if (status) {
        RCLCPP_FATAL(this->get_logger(), "{{ model.name }}acados_create() failed with status %d.", status);
        rclcpp::shutdown();
    }

    sim_config_ = {{ model.name }}_acados_get_sim_config(sim_capsule_);
    sim_dims_ = {{ model.name }}_acados_get_sim_dims(sim_capsule_);
    sim_in_ = {{ model.name }}_acados_get_sim_in(sim_capsule_);
    sim_out_ = {{ model.name }}_acados_get_sim_out(sim_capsule_);
    sim_opts_ = {{ model.name }}_acados_get_sim_opts(sim_capsule_);

    RCLCPP_INFO(this->get_logger(), "acados solver initialized successfully.");
}

void {{ ClassName }}::integration_step() {
    std::array<double, {{ model.name | upper }}_NU> u{};

    {
        std::scoped_lock lock(data_mutex_);
        u = current_u_;
    }

    // Update solver
    this->set_u(u.data());

    // Integrate
    int status = this->sim_solve();
    this->sim_status_behaviour(status);
}

void {{ ClassName }}::sim_status_behaviour(int status) {
    this->get_next_state(xn_.data());
    this->publish_state(xn_, status);
}


// --- ROS Callbacks ---
void {{ ClassName }}::control_callback(const {{ ros_opts.package_name }}_interface::msg::ControlInput::SharedPtr msg) {
    std::scoped_lock lock(data_mutex_);
    std::copy_n(msg->u.begin(), {{ model.name | upper }}_NU, current_u_.begin());
}


// --- ROS Publisher ---
void {{ ClassName }}::publish_state(const std::array<double, {{ model.name | upper }}_NX>& xn, int status) {
    auto state_msg = std::make_unique<{{ ros_opts.package_name }}_interface::msg::State>();
    state_msg->header.stamp = this->get_clock()->now();
    state_msg->header.frame_id = "";
    state_msg->status = status;
    std::copy_n(xn.begin(), {{ model.name | upper }}_NX, state_msg->x.begin());
    state_pub_->publish(std::move(state_msg));
}


// --- Parameter Handling Methods ---
void {{ ClassName }}::setup_parameter_handlers() {
    parameter_handlers_["solver_options.Tsim"] =
        [this](const rclcpp::Parameter& p, rcl_interfaces::msg::SetParametersResult& res) {
            this->set_integration_period(p.as_double());
            // Restart timer with the new period
            if (!this->is_running()) {
                // Assumes that the node is not yet ready and that the timer will be started later
                return;
            }
            try {
                this->start_integration_timer(this->config_.ts);
            } catch (const std::exception& e) {
                res.reason = "Failed to start integration timer, while setting parameter '" + p.get_name() + "': " + e.what();
                res.successful = false;
            }
        };
}

void {{ ClassName }}::declare_parameters() {
    this->declare_parameter("solver_options.Tsim", {{ solver_options.Tsim }});
}

void {{ ClassName }}::load_parameters() {
    this->get_parameter("solver_options.Tsim", config_.ts);
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

// template <size_t N>
// void {{ ClassName }}::get_and_check_array_param(
//     const std::string& param_name,
//     std::array<double, N>& destination
// ) {
//     auto param_value = this->get_parameter(param_name).as_double_array();

//     if (param_value.size() != N) {
//         RCLCPP_ERROR(this->get_logger(), "Parameter '%s' has the wrong size. Expected: %ld, got: %ld",
//                      param_name.c_str(), N, param_value.size());
//         return;
//     }
//     std::copy_n(param_value.begin(), N, destination.begin());
// }

// template <size_t N>
// void {{ ClassName }}::update_param_array(
//     const rclcpp::Parameter& param,
//     std::array<double, N>& destination_array,
//     rcl_interfaces::msg::SetParametersResult& result
// ) {
//     auto values = param.as_double_array();

//     if (values.size() != N) {
//         result.successful = false;
//         result.reason = "Parameter '" + param.get_name() + "' has size " +
//                         std::to_string(values.size()) + ", but expected is " + std::to_string(N) + ".";
//         return;
//     }

//     std::copy_n(values.begin(), N, destination_array.begin());
// }


// --- Helpers ---
void {{ ClassName }}::set_integration_period(double period_seconds) {
    if (config_.ts == period_seconds) {
        // Nothing to do
        return;
    }
    RCLCPP_INFO_STREAM(
        this->get_logger(), "update integration period 'Ts' = " << period_seconds << "s");
    this->config_.ts = period_seconds;
    // Check period validity
    if (this->config_.ts <= 0.0) {
        this->config_.ts = 0.02;
        RCLCPP_WARN(this->get_logger(),
            "Integration period must be positive, defaulting to 0.02s.");
    }
    // Update sim solver integration time
    this->set_t_sim(this->config_.ts);
}

void {{ ClassName }}::start_integration_timer(double period_seconds) {
    this->set_integration_period(period_seconds);
    // if timer already exists, restart with new period
    if (this->is_running()) {
        RCLCPP_WARN(this->get_logger(), "Integration timer already running, restarting...");
        integration_timer_->cancel();
    }
    RCLCPP_INFO_STREAM(this->get_logger(),
        "Starting integration loop with period " << period_seconds << "s.");
    // create timer
    auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::duration<double>(period_seconds));
    integration_timer_ = this->create_wall_timer(
        period,
        std::bind(&{{ ClassName }}::integration_step, this));
}


// --- Acados Helpers ---
int {{ ClassName }}::sim_solve() {
    int status = {{ model.name }}_acados_sim_solve(sim_capsule_);
    if (status != ACADOS_SUCCESS) {
        RCLCPP_ERROR(this->get_logger(), "Simulation Solver failed with status: %d", status);
    }
    return status;
}

void {{ ClassName }}::get_next_state(double* xn) {
    sim_out_get(sim_config_, sim_dims_, sim_out_, "xn", xn);
}

void {{ ClassName }}::set_u(double* u) {
    sim_in_set(sim_config_, sim_dims_, sim_in_, "u", u);
}

void {{ ClassName }}::set_t_sim(double T_sim) {
    sim_in_set(sim_config_, sim_dims_, sim_in_, "T", &T_sim);
}
} // namespace {{ ros_opts.package_name }}


// --- Main ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);
    auto node = std::make_shared<{{ ros_opts.package_name }}::{{ ClassName }}>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
