
#include "acados_c/options.h"

#include "acados_cpp/ocp_qp/options.hpp"

namespace acados {

using std::map;
using std::string;

void flatten(map<string, option_t *>& input, map<string, option_t *>& output) {
    for (auto opt : input) {
        if (opt.second->nested()) {
            for (auto nested_opt : *opt.second) {
                input.erase(opt.first);
                input[opt.first + "." + nested_opt.first] = nested_opt.second;
                flatten(input, output);
            }
        } else {
            output[opt.first] = opt.second;
        }
    }
}

void process_options(string solver_name, map<string, option_t *>& options, void *args) {
    map<string, option_t *> solver_options;
    auto nested_options = std::make_unique<option<map<string, option_t *>>>(options);
    solver_options[solver_name] = nested_options.get();

    auto flattened_options = map<string, option_t *>();
    flatten(solver_options, flattened_options);

    for (auto opt : flattened_options) {
        string option_name = opt.first;
        option_t *opt_p = opt.second;
        bool found = set_option_int(args, option_name.c_str(), std::to_int(opt_p));
        found |= set_option_double(args, option_name.c_str(), std::to_double(opt_p));
        if (!found)
            throw std::invalid_argument("Option '" + option_name + "' is not a valid option.");
    }
}

}  // namespace acados
