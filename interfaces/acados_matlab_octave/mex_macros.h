/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
 */

#include "mex.h"


#define MEX_DIM_CHECK_VEC(fun_name, field, matlab_size, acados_size) {\
    if (acados_size != matlab_size)\
    {\
        sprintf(buffer, "%s: error setting %s, wrong dimension, got %d, need %d", fun_name, field, matlab_size, acados_size);\
        mexErrMsgTxt(buffer);\
    }\
}

#define MEX_DIM_CHECK_VEC_STAGE(fun_name, field, stage, matlab_size, acados_size) {\
    if (acados_size != matlab_size)\
    {\
        sprintf(buffer, "%s: error setting %s for stage %d, wrong dimension, got %d, need %d", fun_name, field, stage, matlab_size, acados_size);\
        mexErrMsgTxt(buffer);\
    }\
}

#define MEX_DIM_CHECK_VEC_TWO(fun_name, field, matlab_size, acados_size1, acados_size2) {\
    if ((acados_size1 != matlab_size) && (acados_size2 != matlab_size) )\
    {\
        sprintf(buffer, "%s: error setting %s, wrong dimension, got %d, need %d or %d",\
                fun_name, field, matlab_size, acados_size1, acados_size2);\
        mexErrMsgTxt(buffer);\
    }\
}

#define MEX_DIM_CHECK_MAT(fun_name, field, matlab_nrow, matlab_ncol, acados_nrow, acados_ncol) {\
    if ((acados_nrow != matlab_nrow) || (acados_ncol != matlab_ncol)  )\
    {\
        sprintf(buffer, "%s: error setting %s, wrong dimension, got %d x %d, need %d x %d",\
               fun_name, field, matlab_nrow, matlab_ncol, acados_nrow, acados_ncol);\
        mexErrMsgTxt(buffer);\
    }\
}

#define MEX_NONBINARY_MAT(fun_name, field) {\
    sprintf(buffer, "%s: error setting %s, matrix should be binary!", fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_MULTIPLE_ONES_IN_ROW(fun_name, field) {\
    sprintf(buffer, "%s: error setting %s, matrix cannot contain multiple ones in row!",\
            fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_SETTER_NO_STAGE_SUPPORT(fun_name, field) {\
    sprintf(buffer, "%s setting %s for specific stage not supported", fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_SETTER_NO_ALL_STAGES_SUPPORT(fun_name, field) {\
    sprintf(buffer, "%s setting %s for all stages is not supported", fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_NOT_SUPPORTED(fun_name, field) {\
    sprintf(buffer, "%s: field %s not supported", fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_NOT_SUPPORTED_SUGGEST(fun_name, field, suggestions) {\
    sprintf(buffer, "%s: field %s not supported, supported values are:\n%s\n", fun_name, field, suggestions);\
    mexErrMsgTxt(buffer);\
}


#define MEX_FIELD_NOT_SUPPORTED_FOR_COST_STAGE(fun_name, field, cost, ii) {\
    sprintf(buffer, "%s: field %s not supported for cost type %d at stage %d", fun_name, field, cost, ii);\
    mexErrMsgTxt(buffer);\
}

#define MEX_MISSING_ARGUMENT(fun_name, field) {\
    sprintf(buffer, "%s: field %s not provided, is mandatory!", fun_name, field);\
    mexErrMsgTxt(buffer);\
}

#define MEX_MISSING_ARGUMENT_NOTE(fun_name, field, note) {\
    sprintf(buffer, "%s: field %s not provided, is mandatory!\nNote: %s", fun_name, field, note);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_VALUE_NOT_SUPPORTED_SUGGEST(fun_name, field, value, suggestions) {\
    sprintf(buffer, "%s: field %s does not support %s, supported values are:\n%s\n",\
            fun_name, field, value, suggestions);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_NOT_SUPPORTED_GIVEN(fun_name, field, value, given, suggestions) {\
    sprintf(buffer, "%s: field %s does not support %s given %s, supported values are:\n%s\n",\
            fun_name, field, value, given, suggestions);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_ONLY_SUPPORTED_FOR_SOLVER(fun_name, field, solver) {\
    sprintf(buffer, "%s: field %s only supported for %s\n",\
            fun_name, field, solver);\
    mexErrMsgTxt(buffer);\
}

#define MEX_FIELD_NOT_SUPPORTED_FOR_SOLVER(fun_name, field, solver) {\
    sprintf(buffer, "%s: field %s not supported for %s\n",\
            fun_name, field, solver);\
    mexErrMsgTxt(buffer);\
}

#define MEX_MISSING_ARGUMENT_MODULE(fun_name, field, module) {\
    sprintf(buffer, "%s field: %s not provided, is mandatory for module %s!", fun_name, field, module);\
    mexErrMsgTxt(buffer);\
}


#define MEX_STR_TO_BOOL(fun_name, str, bool_ptr, name) {\
    if (!strcmp(str, "true"))\
        *bool_ptr = true;\
    else if (!strcmp(str, "false"))\
        *bool_ptr = false;\
    else \
    {\
        sprintf(buffer,\
           "%s error casting field %s to bool, contains %s! possible values \'true\', \'false\'",\
            fun_name, name, str);\
        mexErrMsgTxt(buffer);\
    }\
}


#define MEX_CHECK_DIAGONALITY(fun_name, dim, mat, name){\
    for (int ii=0; ii<dim; ii++)\
    {\
        for (int jj=0; jj<dim; jj++)\
        {\
            if ((jj!=ii) && mat[ii*dim+jj]!= 0.0)\
            {\
                sprintf(buffer,"%s: setting %s, only diagonal matrices supported,\
                            got nonzero offdiagonal elements", fun_name, name);\
                mexErrMsgTxt(buffer);\
            }\
        }\
    }\
}



// macro to string
#define STR(x) STR_AGAIN(x)
#define STR_AGAIN(x) #x

// glue macros
#define GLUE2(x,y) GLUE2_AGAIN(x,y)
#define GLUE2_AGAIN(x,y) x##y
