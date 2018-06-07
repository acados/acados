
%{

#include <string>
#include <vector>

#include "acados_cpp/ocp_qp/ocp_qp.hpp"
#include "acados_cpp/ocp_qp/ocp_qp_solution.hpp"
#include "acados_cpp/options.hpp"

%}

// options.hpp
%rename("$ignore", %$isconstructor) option_t;
%ignore flatten;
%ignore process_options;
%ignore to_double;
%ignore to_int;
%ignore to_map;
%ignore to_string;
%include "acados_cpp/options.hpp"

// ocp_qp_solution.hpp
%rename("$ignore", %$isconstructor) ocp_qp_solution;
%include "acados_cpp/ocp_qp/ocp_qp_solution.hpp"

// ocp_qp.hpp
%rename("$ignore", %$isconstructor) ocp_qp;
%ignore get_field;
%ignore fields;
%include "acados_cpp/ocp_qp/ocp_qp.hpp"

// unignore functions extended below
%rename("%s") ocp_qp;
%rename("%s") get_field;
%rename("%s") fields;

%extend acados::ocp_qp {

    ocp_qp(int N=10, int nx=2, int nu=1) {
        return new acados::ocp_qp(N, nx, nu);
    }

    LangObject *get_field(std::string field) {
        std::vector<std::vector<double>> tmp = $self->get_field(field);
        std::vector<LangObject *> result;
        for (int i = 0; i < tmp.size(); ++i)
            result.push_back(new_matrix($self->shape_of_field(field, i), tmp.at(i).data()));
        return swig::from(result);
    }

    std::vector<std::string> fields() {
        return $self->fields;
    }

}
