#include <vector>

#include "acados_c/ocp_qp.h"

class OcpQp {
    public:
        OcpQp(int N,
              std::vector<int> stateDimensions,
              std::vector<int> controlDimensions,
              std::vector<int> boundDimensions,
              std::vector<int> constraintDimensions);

        const int N;

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

    private:
        void copyDimensions(std::vector<int> dimensions, int *dimension_ptr);

        ocp_qp_dims *dimensions;
        ocp_qp_in *qp;
};
