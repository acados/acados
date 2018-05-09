
#include "acados_cpp/ocp.hpp"

#include <algorithm>

#include "acados/utils/types.h"
#include "acados_cpp/ocp_bounds.hpp"

namespace acados
{
using std::map;
using std::string;
using std::vector;

void ocp::squeeze_dimensions(map<string, vector<vector<double>>> bounds)
{
    // Get bound indices for all stages
    map<string, vector<vector<int>>> idxb_new;
    idxb_new["x"] = calculate_all_idxb(bounds.at("lbx"), bounds.at("ubx"));
    idxb_new["u"] = calculate_all_idxb(bounds.at("lbu"), bounds.at("ubu"));

    // Calculate new dimensions
    map<string, vector<int>> nb{{"x", std::vector<int>(num_stages() + 1)},
                                {"u", std::vector<int>(num_stages() + 1)}};
    std::transform(idxb_new["x"].begin(), idxb_new["x"].end(), nb["x"].begin(),
                   [](const vector<int>& elem) { return elem.size(); });
    std::transform(idxb_new["u"].begin(), idxb_new["u"].end(), nb["u"].begin(),
                   [](const vector<int>& elem) { return elem.size(); });

    change_bound_dimensions(nb["x"], nb["u"]);

    for (string bound : {"x", "u"})
        for (int stage = 0; stage <= num_stages(); ++stage)
            set_bound_indices(bound, stage, idxb_new.at(bound).at(stage));

    needs_initializing(true);
}

void ocp::fill_bounds(map<string, vector<vector<double>>> bounds)
{
    for (auto bound : bounds)
    {
        for (int stage = 0; stage <= num_stages(); ++stage)
        {
            auto idxb_stage = get_bound_indices(bound.first, stage);
            auto stage_bound = bound.second.at(stage);

            vector<double> new_bound(idxb_stage.size());

            if (bound.first == "lbx" || bound.first == "lbu")
            {
                std::fill(new_bound.begin(), new_bound.end(), ACADOS_NEG_INFTY);
            }
            else if (bound.first == "ubx" || bound.first == "ubu")
            {
                std::fill(new_bound.begin(), new_bound.end(), ACADOS_POS_INFTY);
            }

            copy_at(new_bound, stage_bound, idxb_stage);
            set_bound(bound.first, stage, new_bound);
        }
    }
}

}  // namespace acados
