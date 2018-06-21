
#include "acados_cpp/ocp_nlp/casadi_module.hpp"

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <direct.h>
#include <windows.h>
#else
#include <dlfcn.h>
#include <unistd.h>
#endif

#include <utility>

#include "acados_cpp/ocp_nlp/dynamic_loading.hpp"

#include "casadi/mem.h"

namespace acados
{

casadi_module::casadi_module() : function_(), external_function_(), handle_(nullptr, &free_handle)
{

}

casadi_module::casadi_module(const casadi::Function& function, std::string output_folder = "tmp") :
    function_(function), external_function_(), handle_(nullptr, &free_handle)
{
    load_functions(output_folder);
}

casadi_module::casadi_module(casadi_module&& other) :
    function_(),
    external_function_(),
    handle_(nullptr, &free_handle)
{
    *this = std::move(other);
}

casadi_module& casadi_module::operator=(casadi_module&& other)
{
    function_ = other.function_;
    std::swap(handle_, other.handle_);

    // copy external_function_ struct
    external_function_ = other.external_function_;

    // copy double pointers
    for (int i = 0; i < external_function_.args_num; ++i)
        external_function_.args[i] = other.external_function_.args[i];

    for (int i = 0; i < external_function_.res_num; ++i)
        external_function_.res[i] = other.external_function_.res[i];

    other.external_function_.ptr_ext_mem = nullptr;

    generated_header_ = "";
    std::swap(generated_header_, other.generated_header_);

    return *this;
}

casadi_module::~casadi_module()
{
    external_function_casadi_free(&external_function_);
}

const external_function_casadi *casadi_module::as_external_function() const
{
    return &external_function_;
}

void casadi_module::load_functions(std::string output_folder)
{
    generate(output_folder);

    handle_.reset(compile_and_load_library(output_folder, function_.name()));

    external_function_.casadi_fun = reinterpret_cast<casadi_eval_t>(
        load_function(function_.name(), handle_.get()));
    external_function_.casadi_n_in = reinterpret_cast<casadi_getint_t>(
        load_function(function_.name() + "_n_in", handle_.get()));
    external_function_.casadi_n_out = reinterpret_cast<casadi_getint_t>(
        load_function(function_.name() + "_n_out", handle_.get()));
    external_function_.casadi_sparsity_in = reinterpret_cast<casadi_sparsity_t>(
        load_function(function_.name() + "_sparsity_in", handle_.get()));
    external_function_.casadi_sparsity_out = reinterpret_cast<casadi_sparsity_t>(
        load_function(function_.name() + "_sparsity_out", handle_.get()));
    external_function_.casadi_work = reinterpret_cast<casadi_work_t>(
        load_function(function_.name() + "_work", handle_.get()));

    external_function_casadi_create(&external_function_);
}

void casadi_module::generate(std::string output_folder)
{
    casadi::Dict options = {{"with_header", true}, {"with_export", false}, {"casadi_int", "int"}};

    create_directory(output_folder);

    // Hacky workaround because of bug in gcc-4.9

    chdir(output_folder.c_str());

    generated_header_ = output_folder + "/" + function_.generate(options);

    generated_header_.back() = 'h';

    chdir("..");
}

std::string casadi_module::path_to_header()
{
    return generated_header_;
}

std::string casadi_module::name()
{
    return function_.name();
}

}  // namespace acados
