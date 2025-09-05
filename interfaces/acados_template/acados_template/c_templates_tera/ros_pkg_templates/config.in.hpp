#ifndef {{ ros_opts.package_name | upper }}_CONFIG_H
#define {{ ros_opts.package_name | upper }}_CONFIG_H

#include <array>
#include <vector>
#include <string>
#include "acados_solver_{{ model.name }}.h"


namespace {{ ros_opts.package_name }}
{
{%- set ClassName = ros_opts.node_name | replace(from="_", to=" ") | title | replace(from=" ", to="") %}

struct {{ ClassName }}SolverOptions {
    std::string nlp_solver_type{ "{{ solver_options.nlp_solver_type }}" };
    double Tsim{ {{ solver_options.Tsim }} };
};

struct {{ ClassName }}Constraints {
    {%- if dims.nbx > 0 or dims.nbx_e > 0 %}
    // States Bounds
    {%- endif %}
    {%- if dims.nbx > 0 %}
    std::array<double, {{ model.name | upper }}_NBX> lbx{};
    std::array<double, {{ model.name | upper }}_NBX> ubx{};
    {%- endif %}
    {%- if dims.nbx_e > 0 %}
    std::array<double, {{ model.name | upper }}_NBXN> lbx_e{};
    std::array<double, {{ model.name | upper }}_NBXN> ubx_e{};
    {%- endif %}
    {%- if dims.nbu > 0 %}
    
    // Input Bounds
    std::array<double, {{ model.name | upper }}_NBU> lbu{};
    std::array<double, {{ model.name | upper }}_NBU> ubu{};
    {%- endif %}
    {%- if dims.nh_0 > 0 or dims.nh > 0 and dims.nh_e > 0 %}

    // Nonlinear Constraints Bounds (h)
    {%- endif %}
    {%- if dims.nh_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NH0> lh_0{};
    std::array<double, {{ model.name | upper }}_NH0> uh_0{};
    {%- endif %}
    {%- if dims.nh > 0 %}
    std::array<double, {{ model.name | upper }}_NH> lh{};
    std::array<double, {{ model.name | upper }}_NH> uh{};
    {%- endif %}
    {%- if dims.nh_e > 0 %}
    std::array<double, {{ model.name | upper }}_NHN> lh_e{};
    std::array<double, {{ model.name | upper }}_NHN> uh_e{};
    {%- endif %}
    {%- if dims.nphi_0 > 0 or dims.nphi > 0 and dims.nphi_e > 0 %}

    // Nonlinear Phase Constraints Bounds (phi)
    {%- endif %}
    {%- if dims.nphi_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NPHI0> lphi_0{};
    std::array<double, {{ model.name | upper }}_NPHI0> uphi_0{};
    {%- endif %}
    {%- if dims.nphi > 0 %}
    std::array<double, {{ model.name | upper }}_NPHI> lphi{};
    std::array<double, {{ model.name | upper }}_NPHI> uphi{};
    {%- endif %}
    {%- if dims.nphi_e > 0 %}
    std::array<double, {{ model.name | upper }}_NPHIN> lphi_e{};
    std::array<double, {{ model.name | upper }}_NPHIN> uphi_e{};
    {%- endif %}
    {%- if dims.ng > 0 or dims.ng_e > 0 %}

    // General Polytopic Inequalities Bounds (g)
    {%- endif %}
    {%- if dims.ng > 0 %}
    std::array<double, {{ model.name | upper }}_NG> lg{};
    std::array<double, {{ model.name | upper }}_NG> ug{};
    {%- endif %}
    {%- if dims.ng_e > 0 %}
    std::array<double, {{ model.name | upper }}_NGN> lg_e{};
    std::array<double, {{ model.name | upper }}_NGN> ug_e{};
    {%- endif %}
    {%- if dims.nsbx > 0 or dims.nsbx_e > 0 %}

    // State Slack Bounds
    {%- endif %}
    {%- if dims.nsbx > 0 %}
    std::array<double, {{ model.name | upper }}_NSBX> lsbx{};
    std::array<double, {{ model.name | upper }}_NSBX> usbx{};
    {%- endif %}
    {%- if dims.nsbx_e > 0 %}
    std::array<double, {{ model.name | upper }}_NSBXN> lsbx_e{};
    std::array<double, {{ model.name | upper }}_NSBXN> usbx_e{};
    {%- endif %}

    {%- if dims.nsbu > 0 %}
    // Input Slack Bounds
    std::array<double, {{ model.name | upper }}_NSBU> lsbu{};
    std::array<double, {{ model.name | upper }}_NSBU> usbu{};
    {%- endif %}
    {%- if dims.nsh_0 > 0 or dims.nsh > 0 or dims.nsh_e > 0 %}

    // Nonlinear Slack Bounds (sh)
    {%- endif %}
    {%- if dims.nsh_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NSH0> lsh_0{};
    std::array<double, {{ model.name | upper }}_NSH0> ush_0{};
    {%- endif %}
    {%- if dims.nsh > 0 %}
    std::array<double, {{ model.name | upper }}_NSH> lsh{};
    std::array<double, {{ model.name | upper }}_NSH> ush{};
    {%- endif %}
    {%- if dims.nsh_e > 0 %}
    std::array<double, {{ model.name | upper }}_NSHN> lsh_e{};
    std::array<double, {{ model.name | upper }}_NSHN> ush_e{};
    {%- endif %}
    {%- if dims.nsphi_0 > 0 or dims.nsphi > 0 or dims.nsphi_e > 0 %}

    // Nonlinear Phase Slack Bounds (sphi)
    {%- endif %}
    {%- if dims.nsphi_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NSPHI0> lsphi_0{};
    std::array<double, {{ model.name | upper }}_NSPHI0> usphi_0{};
    {%- endif %}
    {%- if dims.nsphi > 0 %}
    std::array<double, {{ model.name | upper }}_NSPHI> lsphi{};
    std::array<double, {{ model.name | upper }}_NSPHI> usphi{};
    {%- endif %}
    {%- if dims.nsphi_e > 0 %}
    std::array<double, {{ model.name | upper }}_NSPHIN> lsphi_e{};
    std::array<double, {{ model.name | upper }}_NSPHIN> usphi_e{};
    {%- endif %}
    {%- if dims.nsg > 0 or dims.nsg_e > 0 %}

    // General Polytopic Slack Bounds (sg)
    {%- endif %}
    {%- if dims.nsg > 0 %}
    std::array<double, {{ model.name | upper }}_NSG> lsg{};
    std::array<double, {{ model.name | upper }}_NSG> usg{};
    {%- endif %}
    {%- if dims.nsg_e > 0 %}
    std::array<double, {{ model.name | upper }}_NSGN> lsg_e{};
    std::array<double, {{ model.name | upper }}_NSGN> usg_e{};
    {%- endif %}
};

struct {{ ClassName }}Weights {
    {%- if dims.ny_0 > 0 %}
    std::array<double, {{ model.name | upper }}_NY0> W_0{};
    {%- endif %}
    {%- if dims.ny > 0 %}
    std::array<double, {{ model.name | upper }}_NY> W{};
    {%- endif %}
    {%- if dims.ny_e > 0 %}
    std::array<double, {{ model.name | upper }}_NYN> W_e{};
    {%- endif %}
};

struct {{ ClassName }}Slacks {
    {%- if dims.ns_0 > 0 %}
    // Initial Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_0')) %}
    std::array<double, {{ model.name | upper }}_NS0> {{ field }}{};
    {%- endif %}
    {%- endfor %}
    
    {%- endif %}
    {%- if dims.ns > 0 %}
    // Stage Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is not ending_with('_0')) and (field is not ending_with('_e')) %}
    std::array<double, {{ model.name | upper }}_NS> {{ field }}{};
    {%- endif %}
    {%- endfor %}
    {%- endif %}
    {%- if dims.ns_e > 0 %}
    
    // Terminal Slacks
    {%- for field, param in cost %}
    {%- set field_l = field | lower %}
    {%- if param and (field_l is starting_with('z')) and (field is ending_with('_e')) %}
    std::array<double, {{ model.name | upper }}_NSN> {{ field }}{};
    {%- endif %}
    {%- endfor %}
    {%- endif %}
};

struct {{ ClassName }}Config {
    {{ ClassName }}SolverOptions    solver_options{};
    {{ ClassName }}Constraints      constraints{};
    {{ ClassName }}Weights          weights{};
    {{ ClassName }}Slacks           slacks{};
    std::array<double, {{ model.name | upper }}_NP>  parameter_values{};
};

} // namespace {{ ros_opts.package_name }}

#endif // {{ ros_opts.package_name | upper }}_CONFIG_H