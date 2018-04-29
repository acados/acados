
#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"

#include <string>

#include "casadi/casadi.hpp"
#include "casadi/mem.h"

namespace acados {

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
std::string compiler {"mex"};
#else
#include <dlfcn.h>
std::string compiler {"cc"};
#endif
int global_library_counter = 0;

static void *load_function(std::string function_name, void *handle) {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
        return GetProcAddress((HMODULE) handle, function_name.c_str());
    #else
        return dlsym(handle, function_name.c_str());
    #endif
}

void ocp_nlp::set_model(const casadi::Function& model, std::map<std::string, option_t *> options) {

    casadi::Dict opts;
    opts["with_header"] = true;
    opts["with_export"] = false;
    model.generate(model.name() + ".c", opts);

    void *handle = compile_and_load(model.name());

    

};

static std::string load_error_message() {
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)

    // Retrieve the system error message for the last-error code
    LPVOID lpMsgBuf;
    DWORD dw = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER |
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL);

    return std::string((LPTSTR) lpMsgBuf);

    #else

    return std::string(dlerror());

    #endif

}

void *compile_and_load(std::string name) {
    void *handle;
    std::string library_name = name + std::to_string(global_library_counter++) + std::string(".so");
    std::string path_to_library = std::string("./") + library_name;
    char command[MAX_STR_LEN];
    snprintf(command, sizeof(command), "%s -fPIC -shared -g %s.c -o %s", compiler.c_str(), name.c_str(),
        library_name.c_str());
    int compilation_failed = system(command);
    if (compilation_failed)
        throw std::runtime_error("Something went wrong when compiling the model.");
    #if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    handle = LoadLibrary(path_to_library.c_str());
    #else
    handle = dlopen(path_to_library.c_str(), RTLD_LAZY);
    #endif
    if (handle == NULL)
        throw std::runtime_error("Loading of " + path_to_library + " failed. Error message: "
                                 + load_error_message());
    return handle;
}

}  // namespace acados
