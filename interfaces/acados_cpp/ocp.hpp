
#ifndef INTERFACES_ACADOS_CPP_OCP_NLP_OCP_H_
#define INTERFACES_ACADOS_CPP_OCP_NLP_OCP_H_

#include <map>
#include <string>
#include <vector>

namespace acados {

class ocp {

 protected:
    
    virtual void squeeze_dimensions(std::map<std::string, std::vector<std::vector<double>>>) = 0;

    virtual void change_bound_dimensions(std::vector<int> nbx, std::vector<int> nbu) = 0;

    virtual bool needs_initializing() = 0;

    virtual void needs_initializing(bool) = 0;

    virtual void set_bounds_indices(std::string, int, std::vector<int>) = 0;

    virtual int num_stages() = 0;

    virtual ~ocp() {};

};

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OCP_NLP_OCP_H_
