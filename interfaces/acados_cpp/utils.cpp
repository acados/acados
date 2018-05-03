
#include "acados_cpp/utils.hpp"

namespace acados {

std::string to_string(std::pair<int, int> p) {
    return "( " + std::to_string(p.first) + ", " + std::to_string(p.second) + " )";
}

template<typename T>
std::string to_string(std::vector<T> v) {
    std::string result_string = " vector of length " + std::to_string(v.size()) + ": [\n ";
    for(auto it : v) {
        result_string += std::to_string(it) + ", ";
    }
    return result_string + "]\n";
}

bool match(std::pair<int, int> dims, int nb_elems) {
    int nb_expected_elems = dims.first * dims.second;
    if (nb_expected_elems == 0 || nb_expected_elems == nb_elems)
        return true;
    return false;
}

template<typename T>
const T& clamp(const T& lo, const T& hi, const T& val) {
    if (val < lo)
        return lo;
    else if (val > hi)
        return hi;
    
    return val;
}

std::string load_error_message() {
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

}  // namespace acados
