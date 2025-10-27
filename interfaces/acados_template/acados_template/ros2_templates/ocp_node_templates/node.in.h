#ifndef {{ ros_opts.node_name | upper }}_H
#define {{ ros_opts.node_name | upper }}_H

#include <rclcpp/rclcpp.hpp>
#include <rcutils/logging.h>
#include <mutex>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

// ROS2 message includes
#include "{{ ros_opts.package_name }}_interface/msg/state.hpp"
#include "{{ ros_opts.package_name }}_interface/msg/control.hpp"
#include "{{ ros_opts.package_name }}_interface/msg/references.hpp"
{%- if ros_opts.publish_control_sequence %}
#include "{{ ros_opts.package_name }}_interface/msg/control_sequence.hpp"
{%- endif %}
{%- if dims.np > 0 %}
#include "{{ ros_opts.package_name }}_interface/msg/parameters.hpp"
{%- endif %}
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
{%- set use_multithreading = ros_opts.threads is defined and ros_opts.threads > 1 %}
class {{ ClassName }} : public rclcpp::Node {
private:
    // --- ROS Subscriptions ---
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::State>::SharedPtr state_sub_;
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::References>::SharedPtr references_sub_;
    {%- if dims.np > 0 %}
    rclcpp::Subscription<{{ ros_opts.package_name }}_interface::msg::Parameters>::SharedPtr parameters_sub_;
    {%- endif %}

    // --- ROS Publishers ---
    rclcpp::Publisher<{{ ros_opts.package_name }}_interface::msg::Control>::SharedPtr control_pub_;
    {%- if ros_opts.publish_control_sequence %}
    rclcpp::Publisher<{{ ros_opts.package_name }}_interface::msg::ControlSequence>::SharedPtr control_sequence_pub_;
    {%- endif %}

    // --- ROS Params and Timer
    rclcpp::TimerBase::SharedPtr control_timer_;
    OnSetParametersCallbackHandle::SharedPtr param_callback_handle_;
    using ParamHandler = std::function<void(const rclcpp::Parameter&, rcl_interfaces::msg::SetParametersResult&)>;
    std::unordered_map<std::string, ParamHandler> parameter_handlers_;

    // --- Acados Solver ---
    {{ model.name }}_solver_capsule *ocp_capsule_;
    ocp_nlp_in* ocp_nlp_in_;
    ocp_nlp_out* ocp_nlp_out_;
    ocp_nlp_out* ocp_nlp_sens_;
    ocp_nlp_config* ocp_nlp_config_;
    void* ocp_nlp_opts_;
    ocp_nlp_dims* ocp_nlp_dims_;
    {%- if use_multithreading %}

    // --- Multithreading ---
    rclcpp::CallbackGroup::SharedPtr timer_group_;
    rclcpp::CallbackGroup::SharedPtr services_group_;
    std::mutex data_mutex_;
    std::recursive_mutex solver_mutex_;
    {%- endif %}

    // --- Data and States ---
    {{ ClassName }}Config config_;
    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    bool first_solve_{true};
    {%- endif %}
    std::array<double, {{ model.name | upper }}_NU> u0_;
    {%- if ros_opts.publish_control_sequence %}
    std::array<std::array<double, {{ model.name | upper }}_NU>, {{ solver_options.N_horizon }}> u_seq_{};
    {%- endif %}
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
    void solver_status_behaviour(int status);

    // --- ROS Callbacks ---
    void state_callback(const {{ ros_opts.package_name }}_interface::msg::State::SharedPtr msg);
    void references_callback(const {{ ros_opts.package_name }}_interface::msg::References::SharedPtr msg);
    {%- if dims.np > 0 %}
    void parameters_callback(const {{ ros_opts.package_name }}_interface::msg::Parameters::SharedPtr msg);
    {%- endif %}

    // --- ROS Publisher ---
    void publish_control(
        const std::array<double, {{ model.name | upper }}_NU>& u0, 
        int status);
    {%- if ros_opts.publish_control_sequence %}
    void publish_control_sequence(
        const std::array<std::array<double, {{ model.name | upper }}_NU>, {{ solver_options.N_horizon }}>& u_sequence, 
        int status);
    {%- endif %}

    // --- ROS Parameter ---
    void setup_parameter_handlers();
    void declare_parameters();
    void load_parameters();
    void apply_all_parameters_to_solver();
    rcl_interfaces::msg::SetParametersResult on_parameter_update(
        const std::vector<rclcpp::Parameter>& params);

    template <size_t N>
    void get_and_check_array_param(
        const std::string& param_name,
        std::array<double, N>& destination);
    template <size_t N>
    void update_param_array(
        const rclcpp::Parameter& param,
        std::array<double, N>& destination_array,
        rcl_interfaces::msg::SetParametersResult& result);
    template<size_t N>
    void update_constraint(
        const rclcpp::Parameter& param,
        rcl_interfaces::msg::SetParametersResult& result,
        const char* field,
        const std::vector<int>& stages);
    template<size_t N>
    void update_cost(
        const rclcpp::Parameter& param,
        rcl_interfaces::msg::SetParametersResult& result,
        const char* field,
        const std::vector<int>& stages);

    // --- Helpers ---
    void set_period(double period_seconds);
    void start_control_timer(double period_seconds = {{ solver_options.time_steps[0] }});
    bool is_running() const {
        return control_timer_ && !control_timer_->is_canceled();
    }

    // --- Acados Solver ---
    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    int prepare_rti_solve();
    int feedback_rti_solve();
    {%- endif %}
    int ocp_solve();

    // --- Acados Getter ---
    void get_control(double* u, int stage);
    void get_state(double* x, int stage);

    // --- Acados Setter ---
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
    {%- if solver_options.nlp_solver_type == "SQP_RTI" %}
    void warmstart_solver_states(double *x0);
    {%- endif %}

    // --- Helpers ---
    bool check_acados_status(
        const char* field,
        int stage,
        int status);
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.node_name | upper }}_H
