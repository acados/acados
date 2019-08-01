
#include "acados_cpp/options.hpp"

#include <memory>

#include "acados_c/options_interface.h"

namespace acados
{
using std::map;
using std::string;

void flatten(const map<string, option_t *> &input, map<string, option_t *> &output)
{
    for (auto opt : input)
    {
        if (opt.second->nested())
        {
            auto option_name = opt.first;
            for (auto nested_option : opt.second->as_map())
                flatten({{option_name + "." + nested_option.first, nested_option.second}}, output);
        }
        else
        {
            output[opt.first] = opt.second;
        }
    }
}

void process_options(string solver_name, map<string, option_t *> &options, void *args)
{
    map<string, option_t *> solver_options;
    auto nested_options = std::unique_ptr<option<map<string, option_t *>>>(
        new option<map<string, option_t *>>(options));
    solver_options[solver_name] = nested_options.get();

    auto flattened_options = map<string, option_t *>();
    flatten(solver_options, flattened_options);

    for (auto opt : flattened_options)
    {
        string option_name = opt.first;
        option_t *opt_p = opt.second;
        bool found = set_option_int(args, option_name.c_str(), to_int(opt_p));
        found |= set_option_double(args, option_name.c_str(), to_double(opt_p));
        if (!found)
            throw std::invalid_argument("Option '" + option_name + "' is not a valid option.");
    }
}

}  // namespace acados
