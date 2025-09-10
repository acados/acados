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
#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{ model.name }}.h"

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
    rclcpp::TimerBase::SharedPtr control_timer_;
    OnSetParametersCallbackHandle::SharedPtr param_callback_handle_;
    using ParamHandler = std::function<void(const rclcpp::Parameter&, rcl_interfaces::msg::SetParametersResult&)>;
    std::unordered_map<std::string, ParamHandler> parameter_handlers_;

    // --- Acados Solver ---
    {{ model.name }}_solver_capsule *ocp_capsule_;
    ocp_nlp_config* ocp_nlp_config_;
    ocp_nlp_dims* ocp_nlp_dims_;
    ocp_nlp_in* ocp_nlp_in_;
    ocp_nlp_out* ocp_nlp_out_;
    void* ocp_nlp_opts_;

    // --- Data and States ---
    std::mutex data_mutex_;
    {{ ClassName }}Config config_;

    std::array<double, {{ model.name | upper }}_NU> u0_;
    std::array<double, {{ model.name | upper }}_NX> current_x_;

public:
    {{ ClassName }}();
    ~{{ ClassName }}();

private:
    // --- Core Methods ---
    void initialize_simulator();
    void simulation_loop();
    void sim_status_behaviour(int status);

    // --- ROS Callbacks ---
    void control_callback(const {{ ros_opts.package_name }}_interface::msg::ControlInput::SharedPtr msg);

    // --- ROS Publisher ---
    void publish_state(const std::array<double, {{ model.name | upper }}_NX>& x0, int status);

    // --- Parameter Handling Methods ---
    void setup_parameter_handlers();
    void declare_parameters();
    void load_parameters();
    rcl_interfaces::msg::SetParametersResult on_parameter_update(const std::vector<rclcpp::Parameter>& params);

    // --- Helpers ---
    void start_simulation_timer(double period_seconds = 0.02);

    // --- Acados Helpers ---
    int sim_solve();

    void get_input(double* u, int stage);
    void get_state(double* x, int stage);
    void set_x0(double* x0);
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.node_name | upper }}_H