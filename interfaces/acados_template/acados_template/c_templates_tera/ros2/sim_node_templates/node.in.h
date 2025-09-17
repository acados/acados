#ifndef {{ ros_opts.node_name | upper }}_H
#define {{ ros_opts.node_name | upper }}_H

#include <rclcpp/rclcpp.hpp>
#include <mutex>
#include <array>
#include <vector>
#include <unordered_map>

// ROS2 message includes
#include "{{ ros_opts.package_name }}_interface/msg/state.hpp"
#include "{{ ros_opts.package_name }}_interface/msg/control_input.hpp"
#include "std_msgs/msg/header.hpp"

// Acados includes
#include "acados/sim/sim_common.h"
#include "acados_c/sim_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_sim_solver_{{ model.name }}.h"

// Package includes
#include "{{ ros_opts.package_name }}/utils.hpp"
#include "{{ ros_opts.package_name }}/config.hpp"


namespace {{ ros_opts.package_name }}
{

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
class {{ ClassName }} : public rclcpp::Node {
private:
    // --- ROS Subscriptions ---
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::ControlInput>::SharedPtr control_sub_;

    // --- ROS Publishers ---
    rclcpp::Publisher<{{ ros_opts.package_name }}_interface::msg::State>::SharedPtr state_pub_;

    // --- ROS Params and Timer
    rclcpp::TimerBase::SharedPtr integration_timer_;
    OnSetParametersCallbackHandle::SharedPtr param_callback_handle_;
    using ParamHandler = std::function<void(const rclcpp::Parameter&, rcl_interfaces::msg::SetParametersResult&)>;
    std::unordered_map<std::string, ParamHandler> parameter_handlers_;

    // --- Acados Solver ---
    {{ model.name }}_sim_solver_capsule *sim_capsule_;
    sim_config* sim_config_;
    void* sim_dims_;
    sim_in* sim_in_;
    sim_out* sim_out_;
    sim_opts* sim_opts_;

    // --- Data and States ---
    std::mutex data_mutex_;
    {{ ClassName }}Config config_;

    std::array<double, {{ model.name | upper }}_NX> xn_;
    std::array<double, {{ model.name | upper }}_NU> current_u_;

public:
    {{ ClassName }}();
    ~{{ ClassName }}();

private:
    // --- Core Methods ---
    void initialize_simulator();
    void integration_step();
    void sim_status_behaviour(int status);

    // --- ROS Callbacks ---
    void control_callback(const {{ ros_opts.package_name }}_interface::msg::ControlInput::SharedPtr msg);

    // --- ROS Publisher ---
    void publish_state(const std::array<double, {{ model.name | upper }}_NX>& xn, int status);

    // --- Parameter Handling Methods ---
    void setup_parameter_handlers();
    void declare_parameters();
    void load_parameters();
    rcl_interfaces::msg::SetParametersResult on_parameter_update(const std::vector<rclcpp::Parameter>& params);

    // --- Helpers ---
    void set_integration_period(double period_seconds);
    void start_integration_timer(double period_seconds = 0.02);
    bool is_running() const {
        return integration_timer_ && !integration_timer_->is_canceled();
    }

    // --- Acados Helpers ---
    int sim_solve();
    void get_next_state(double* xn);
    void set_u(double* u);
    void set_t_sim(double T_sim);
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.node_name | upper }}_H