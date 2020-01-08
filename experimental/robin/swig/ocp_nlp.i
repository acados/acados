
%{

#include "acados_cpp/ocp_nlp/ocp_nlp.hpp"
#include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"

%}

%rename("$ignore", %$isconstructor) ocp_nlp_solution;
%include "acados_cpp/ocp_nlp/ocp_nlp_solution.hpp"

%include "acados_cpp/ocp_nlp/ocp_nlp.hpp"
