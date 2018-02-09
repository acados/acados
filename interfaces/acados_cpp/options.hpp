
#ifndef ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
#define ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_

#include <map>
#include <string>

namespace acados {

class option_t {
public:
    inline virtual std::string repr() {
        return "";
    }
    
    inline virtual int as_int() {
        return 0;
    }

    inline virtual double as_double() {
        return 0.0;
    }

    inline virtual std::map<std::string, option_t *>& as_map() {
        throw std::runtime_error("option_t is no map");
    }

    inline virtual bool nested() {
        return false;
    }

    inline virtual std::map<std::string, option_t *>::iterator begin() {
        throw std::range_error("option_t is not iterable");
    }

    inline virtual std::map<std::string, option_t *>::iterator end() {
        throw std::range_error("option_t is not iterable");
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

    inline int as_int() override {
        return static_cast<int>(value);
    }

    inline double as_double() override {
        return static_cast<double>(value);
    }

    ~option() override = default;
private:
    T value;
};

template<>
class option<std::map<std::string, option_t *>> : public option_t {
public:
    option(std::map<std::string, option_t *> val) : value(val) {}

    inline bool nested() override {
        return true;
    }

    inline std::map<std::string, option_t *>& as_map() override {
        return value;
    }

    inline std::map<std::string, option_t *>::iterator begin() override {
        return value.begin();
    }

    inline std::map<std::string, option_t *>::iterator end() override {
        return value.end();
    }

private:
    std::map<std::string, option_t *> value;
};

}  // namespace acados

namespace std {

inline string to_string(acados::option_t *opt) {
    return opt->repr();
}

inline int to_int(acados::option_t *opt) {
    return opt->as_int();
}

inline double to_double(acados::option_t *opt) {
    return opt->as_double();
}

inline std::map<std::string, acados::option_t *>& to_map(acados::option_t *opt) {
    return opt->as_map();
}

}  // namespace std

#endif  // ACADOS_INTERFACES_ACADOS_CPP_OPTIONS_HPP_
