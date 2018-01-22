
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "acados_c/ocp_qp.h"

using std::vector;
using std::string;

namespace acados {

class ocp_qp {

public:

    ocp_qp(int N, vector<int> nx, vector<int> nu, vector<int> nbx, vector<int> nbu, vector<int> ng);

    ocp_qp(int N, int nx, int nu, int nbx = 0, int nbu = 0, int ng = 0);

    void update(string field, int stage, vector<double> v);
    vector< vector<double> > extract(string field);

    vector<int> nx();
    vector<int> nu();
    vector<int> nbx();
    vector<int> nbu();
    vector<int> ng();

    friend std::ostream& operator<<(std::ostream& oss, const ocp_qp& qp);

    const int N;

private:
    
    std::pair<int, int> dimensions(string field, int stage);

    void write_dimensions(vector<int> dims, int *ptr);

    void check_range(string field, int stage);
    
    void check_nb_elements(string, int stage, int nb_elems);

    static std::map<string, std::function<void(int, double *, ocp_qp_in *)>> extract_functions;

    ocp_qp_in *qp;
};

}  // namespace acados
