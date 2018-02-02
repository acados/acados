
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_

#include <map>

namespace acados {

class option_t {};

template<typename T>
class option : public option_t {
public:
    option(T val) : value(val) {}
    // option(T *val) : T(*val) {}
private:
    T value;
};

}  // namespace acados

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
