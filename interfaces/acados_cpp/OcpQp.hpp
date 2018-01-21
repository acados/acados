
// #include <utility>
#include <vector>

#include "acados_c/ocp_qp.h"

namespace acados {

class OcpQp {

    public:

        enum CostMatrix {Q, S, R};
        enum CostVector {q, r};
        enum DynamicsMatrix {A, B};
        enum DynamicsVector {b};
        enum Bound {lbx, ubx, lbu, ubu, lg, ug};
        enum ConstraintMatrix {C, D};

        OcpQp(int N,
              std::vector<int> stateDimensions,
              std::vector<int> controlDimensions,
              std::vector<int> stateBoundDimensions,
              std::vector<int> controlBoundDimensions,
              std::vector<int> constraintDimensions);

        OcpQp(int N,
              int stateDimensions,
              int controlDimensions,
              int stateBoundDimensions = 0,
              int controlBoundDimensions = 0,
              int constraintDimensions = 0);

        const int N;

        std::vector<double> getA(int i);
        int numRowsA(int i);
        int numColsA(int i);

        int numRowsB(int i);
        int numColsB(int i);

        int numRowsb(int i);
        int numColsb(int i);

        int numRowsQ(int i);
        int numColsQ(int i);

        int numRowsR(int i);
        int numColsR(int i);

        int numRowsS(int i);
        int numColsS(int i);

        int numRowsq(int i);
        int numColsq(int i);

        int numRowsr(int i);
        int numColsr(int i);

        int numRowslb(int i);
        int numColslb(int i);

        int numRowsub(int i);
        int numColsub(int i);

        int numRowsC(int i);
        int numColsC(int i);

        int numRowsD(int i);
        int numColsD(int i);

        int numRowslg(int i);
        int numColslg(int i);

        int numRowsug(int i);
        int numColsug(int i);

        void setQ(int i, double *Q, int num_rows, int num_cols);
        void setQ(int i, double *Q);
        
        void setS(int i, double *S, int num_rows, int num_cols);
        void setS(int i, double *S);
        
        void setR(int i, double *R, int num_rows, int num_cols);
        void setR(int i, double *R);
        
        void setq(int i, double *q, int num_elems);
        void setq(int i, double *q);

        void setr(int i, double *r, int num_elems);
        void setr(int i, double *r);

        void setA(int i, double *A, int num_rows, int num_cols);
        void setA(int i, double *A);
    
        void setB(int i, double *B, int num_rows, int num_cols);
        void setB(int i, double *B);
        
        void setb(int i, double *b, int num_elems);
        void setb(int i, double *b);

        void setlb(int i, double *lb, int num_elems);
        void setlb(int i, double *lb);

        void setub(int i, double *ub, int num_elems);
        void setub(int i, double *ub);

        void setC(int i, double *C, int num_rows, int num_cols);
        void setC(int i, double *C);
        
        void setD(int i, double *D, int num_rows, int num_cols);
        void setD(int i, double *D);

        void setlg(int i, double *lg, int num_elems);
        void setlg(int i, double *lg);

        void setug(int i, double *ug, int num_elems);
        void setug(int i, double *ug);


        friend std::ostream& operator<<(std::ostream& oss, const OcpQp& qp);

        void copyDimensions(std::vector<int> dimensions, int *dimension_ptr);

        std::pair<int, int> dims(Field elem, int stage);
        std::pair<int, int> offset(Field elem, int stage);

        void update(int stage, Field elem, std::vector<double> v);

        ocp_qp_dims *dimensions;
        ocp_qp_in *qp;
};

}  // namespace acados
