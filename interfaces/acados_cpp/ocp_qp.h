#include <vector>

#include "acados_c/ocp_qp.h"

class OcpQp {
    public:
        OcpQp(int N, int numStates, int numControls, int numBounds = 0, int numConstraints = 0);
        const int N;

    private:
        void setStateDimensionAt(int index, int numStates);
        void setStateDimensions(int numStates);
        void setControlDimensionAt(int index, int numControls);
        void setControlDimensions(int numControls);

        std::vector<int> nx;
        std::vector<int> nu;
        std::vector<int> nb;
        std::vector<int> nc;

        ocp_qp_dims *dimensions;
        ocp_qp_in *qp;
};
