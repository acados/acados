
#include "acados_cpp/ocp_nlp/dynamic_loading.hpp"

#include <cstdlib>
#include <stdexcept>

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace acados {

using std::string;

void *compile_and_load_library(string output_folder, string source_name)
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    std::string dynamic_library_suffix {".dll"};
    std::string compiler {"gcc"};
#else
    std::string dynamic_library_suffix {".so"};
    std::string compiler {"cc"};
#endif
    void *handle;
    std::string library_name = source_name + dynamic_library_suffix;
    std::string path_to_library = "./" + output_folder + "/" + library_name;
    std::string path_to_file = "./" + output_folder + "/" + source_name;
    char command[256];
    snprintf(command, sizeof(command), "%s -fPIC -shared -O3 %s.c -o %s", compiler.c_str(),
             path_to_file.c_str(), path_to_library.c_str());
    int compilation_failed = system(command);
    if (compilation_failed)
        throw std::runtime_error("Something went wrong when compiling the model.");
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    handle = LoadLibrary(path_to_library.c_str());
#else
    handle = dlopen(path_to_library.c_str(), RTLD_LAZY);
#endif
    if (handle == NULL)
        throw std::runtime_error("Loading of " + path_to_library +
                                 " failed. Error message: " + load_error_message());
    return handle;
}

void *load_function(std::string function_name, void *handle)
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    return (void *) GetProcAddress((HMODULE) handle, function_name.c_str());
#else
    return dlsym(handle, function_name.c_str());
#endif
}

void free_handle(void *handle)
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)
    FreeLibrary((HMODULE) handle);
#else
    dlclose(handle);
#endif
}

std::string load_error_message()
{
#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)

    // Retrieve the system error message for the last-error code
    LPVOID lpMsgBuf;
    DWORD dw = GetLastError();

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL, dw, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf, 0, NULL);

    return std::string((LPTSTR) lpMsgBuf);

#else

    return std::string(dlerror());

#endif
}

#if (defined _WIN32 || defined _WIN64 || defined __MINGW32__ || defined __MINGW64__)

void create_directory(std::string name)
{
    try {
        CreateDirectory(name.c_str(), NULL);
    } catch (...) {
        throw;
    }
}

#else

void create_directory(std::string name)
{
    std::string command = "mkdir -p " + name;
    int creation_failed = system(command.c_str());
    if (creation_failed)
        throw std::runtime_error("Could not create folder '" + name + "'.");
}

#endif

}  // namespace acados
