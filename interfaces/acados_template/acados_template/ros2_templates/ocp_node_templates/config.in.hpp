#ifndef {{ ros_opts.package_name | upper }}_CONFIG_H
#define {{ ros_opts.package_name | upper }}_CONFIG_H

#include <array>
#include <vector>
#include <string>


namespace {{ ros_opts.package_name }}
{
{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}

struct {{ ClassName }}Config {
    double ts{ {{ solver_options.Tsim }} };
    int threads{ {{ ros_opts.threads | default(value=1) }} };
    bool verbose{ false };
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.package_name | upper }}_CONFIG_H
