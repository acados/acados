
#include "acados_cpp/ocp.hpp"

#include "acados_cpp/ocp_bounds.hpp"

namespace acados {

using std::map;
using std::string;
using std::vector;

void ocp::squeeze_dimensions(map<string, vector<vector<double>>> bounds) {
    
    // Get bound indices for all stages
    map<string, vector<vector<int>>> idxb_new;
    idxb_new["x"] = calculate_all_idxb(bounds.at("lbx"), bounds.at("ubx"));
    idxb_new["u"] = calculate_all_idxb(bounds.at("lbu"), bounds.at("ubu"));

    // Calculate new dimensions
    map<string, vector<int>> nb {{"x", std::vector<int>(num_stages()+1)}, {"u", std::vector<int>(num_stages()+1)}};
    std::transform(idxb_new["x"].begin(), idxb_new["x"].end(), nb["x"].begin(), [](auto elem){return elem.size();});
    std::transform(idxb_new["u"].begin(), idxb_new["u"].end(), nb["u"].begin(), [](auto elem){return elem.size();});

    change_bound_dimensions(nb["x"], nb["u"]);

    for (string bound : {"x", "u"})
        for (int stage = 0; stage <= num_stages(); ++stage)
            set_bounds_indices(bound, stage, idxb_new.at(bound).at(stage));

    needs_initializing(true);
}

}  // namespace acados
