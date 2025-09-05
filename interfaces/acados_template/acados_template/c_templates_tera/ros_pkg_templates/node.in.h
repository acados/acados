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
#include "{{ ros_opts.package_name }}_interface/msg/references.hpp"
#include "{{ ros_opts.package_name }}_interface/msg/parameters.hpp"
#include "std_msgs/msg/header.hpp"

// Acados includes
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados/ocp_nlp/ocp_nlp_sqp_rti.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "blasfeo_d_aux_ext_dep.h"
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
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::State>::SharedPtr state_sub_;
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::References>::SharedPtr references_sub_;
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::Parameters>::SharedPtr parameters_sub_;

    // --- ROS Publishers ---
    rclcpp::Publisher<{{ ros_opts.package_name }}_interface::msg::ControlInput>::SharedPtr control_input_pub_;

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
    std::array<double, {{ model.name | upper }}_NU> u0_default_;
    std::array<double, {{ model.name | upper }}_NX> current_x_;
    {%- if dims.ny_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NY0> current_yref_0_;
    {%- endif %}
    {%- if dims.ny > 0 %}
    std::array<double, {{ model.name | upper }}_NY> current_yref_;
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    std::array<double, {{ model.name | upper }}_NYN> current_yref_e_;
    {%- endif %}
    {%- if dims.np > 0 %}
    std::array<double, {{ model.name | upper }}_NP> current_p_;
    {%- endif %}

public:
    {{ ClassName }}();
    ~{{ ClassName }}();

private:
    // --- Core Methods ---
    void initialize_solver();
    void control_loop();

    // --- ROS Callbacks ---
    void state_callback(const {{ ros_opts.package_name }}_interface::msg::State::SharedPtr msg);
    void references_callback(const {{ ros_opts.package_name }}_interface::msg::References::SharedPtr msg);
    void parameters_callback(const {{ ros_opts.package_name }}_interface::msg::Parameters::SharedPtr msg);

    // --- ROS Publisher ---
    void publish_input(const std::array<double, {{ model.name | upper }}_NU>& u0);

    // --- Parameter Handling Methods ---
    void setup_parameter_handlers();
    void declare_parameters();
    void load_parameters();
    void log_parameters();
    rcl_interfaces::msg::SetParametersResult on_parameter_update(const std::vector<rclcpp::Parameter>& params);

    // --- Helpers ---
    void start_control_timer(double rate_hz = 50.0);

    // --- Acados Helpers ---
    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    int prepare_rti_solve();
    int feedback_rti_solve();
    {%- endif %}
    int ocp_solve();

    void get_input(double* u, int stage);
    void get_state(double* x, int stage);
    
    void set_x0(double* x0);
    {%- if dims.ny_0 > 0 %}
    void set_yref0(double* yref0);
    {%- endif %}
    {%- if dims.ny > 0 %}
    void set_yref(double* yref, int stage);
    void set_yrefs(double* yref);
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    void set_yref_e(double* yref_e);
    {%- endif %}
    {%- if dims.np > 0 %}
    void set_ocp_parameter(double* p, size_t np, int stage);
    void set_ocp_parameters(double* p, size_t np);
    {%- endif %}

    void set_cost_weights();
    {%- if dims.ns_0 > 0 or dims.ns > 0 or dims.ns_e > 0 %}
    void set_slack_weights();
    {%- endif %}
    void set_constraints();
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.node_name | upper }}_H