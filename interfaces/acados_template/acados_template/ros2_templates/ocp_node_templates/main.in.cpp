#include <rclcpp/version.h>
#include <rcutils/logging.h>

#include "{{ ros_opts.package_name }}/node.h"

{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}
{%- set use_multithreading = ros_opts.threads is defined and ros_opts.threads > 1 %}

// --- Main ---
int main(int argc, char **argv) {
    rclcpp::init(argc, argv);

    // Suppress ROS timer logging
#if RCLCPP_VERSION_MAJOR >= 32
    const char *logger_name = "rclcpp";
#else
    const char *logger_name = "rcl";
#endif
    if (rcutils_logging_set_logger_level(logger_name, RCUTILS_LOG_SEVERITY_WARN) != RCUTILS_RET_OK) {
        RCLCPP_WARN(rclcpp::get_logger(logger_name), "Failed to set logger level for '%s'", logger_name);
    }

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
