
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_

#include <string>

namespace acados {

class option_t {
public:
    inline virtual std::string repr() {
        return "";
    }
    virtual ~option_t() = default;
};

template<typename T>
class option : public option_t {
public:
    option(T val) : value(val) {}
    inline std::string repr() override {
        return std::to_string(value);
    }
    ~option() override = default;
private:
    T value;
};

}  // namespace acados

namespace std {

inline string to_string(acados::option_t *opt) {
    return opt->repr();
}

}  // namespace std

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
