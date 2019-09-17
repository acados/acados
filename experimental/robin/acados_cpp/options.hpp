/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef INTERFACES_ACADOS_CPP_OPTIONS_HPP_
#define INTERFACES_ACADOS_CPP_OPTIONS_HPP_

#include <map>
#include <string>

namespace acados
{
class option_t
{
 public:
    inline virtual std::string repr() { return ""; }

    inline virtual int as_int() { return 0; }

    inline virtual double as_double() { return 0.0; }

    inline virtual std::map<std::string, option_t *> &as_map()
    {
        throw std::runtime_error("option_t is no map");
    }

    inline virtual bool nested() { return false; }

    inline virtual std::map<std::string, option_t *>::iterator begin()
    {
        throw std::range_error("option_t is not iterable");
    }

    inline virtual std::map<std::string, option_t *>::iterator end()
    {
        throw std::range_error("option_t is not iterable");
    }

    virtual ~option_t() = default;
};

void flatten(const std::map<std::string, option_t *> &input,
             std::map<std::string, option_t *> &output);

void process_options(std::string solver_name, std::map<std::string, option_t *> &options,
                     void *args);

template <typename T>
class option : public option_t
{
 public:
    explicit option(T val) : value(val) {}

    inline std::string repr() override { return std::to_string(value); }

    inline int as_int() override { return static_cast<int>(value); }

    inline double as_double() override { return static_cast<double>(value); }

    ~option() override = default;

 private:
    T value;
};

template <>
class option<std::string> : public option_t
{
 public:
    explicit option(std::string val) : value(val) {}

    inline std::string repr() override { return value; }

    inline int as_int() override { return 0; }

    inline double as_double() override { return 0; }

    ~option() override = default;

 private:
    std::string value;
};

template <>
class option<std::map<std::string, option_t *>> : public option_t
{
 public:
    explicit option(std::map<std::string, option_t *> val) : value(val) {}

    inline bool nested() override { return true; }

    inline std::map<std::string, option_t *> &as_map() override { return value; }

    inline std::map<std::string, option_t *>::iterator begin() override { return value.begin(); }

    inline std::map<std::string, option_t *>::iterator end() override { return value.end(); }

 private:
    std::map<std::string, option_t *> value;
};

inline std::string to_string(acados::option_t *opt) { return opt->repr(); }

inline int to_int(acados::option_t *opt) { return opt->as_int(); }

inline double to_double(acados::option_t *opt) { return opt->as_double(); }

inline std::map<std::string, acados::option_t *> &to_map(acados::option_t *opt)
{
    return opt->as_map();
}

}  // namespace acados

#endif  // INTERFACES_ACADOS_CPP_OPTIONS_HPP_
