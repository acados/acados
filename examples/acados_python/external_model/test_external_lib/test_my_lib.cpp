// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
#include "include/external_ode_casadi.h"

int main(){
    // Use CasADi's "external" to load the compiled function
    casadi::Function f = casadi::external(libexternal_ode_casadi_name());
    casadi::DM x(reshape(casadi::DM(std::vector<double>({1, 2})), 2, 1));
    casadi::DM u(reshape(casadi::DM(std::vector<double>({3})), 1, 1));
    // Use like any other CasADi function
    std::vector<casadi::DM> arg = {x, u};
    std::vector<casadi::DM> fRes = f(arg);
    std::cout << "x: "           << x << std::endl;
    std::cout << "u: "           << u << std::endl;
    std::cout << "xDot=f(x,u): " << fRes.at(0) << std::endl;
    return 0;
}
