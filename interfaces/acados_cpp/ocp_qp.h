#include <vector>

#include "acados_c/ocp_qp.h"

class OcpQp {
    public:
        OcpQp(int N,
              std::vector<int> stateDimensions,
              std::vector<int> controlDimensions,
              std::vector<int> boundDimensions,
              std::vector<int> constraintDimensions);

        friend std::ostream& operator<<(std::ostream& oss, const OcpQp& qp);

    private:
        void copyDimensions(std::vector<int> dimensions, int *dimension_ptr);

        ocp_qp_dims *dimensions;
        ocp_qp_in *qp;
};
