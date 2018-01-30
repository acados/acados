/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/utils/external_function.h"



int external_function_dims_calculate_size(int num_inputs, int num_outputs)
{
    int size = sizeof(external_function_dims);

    size += num_inputs * sizeof(int);

    size += num_outputs * sizeof(int);

    return size;
}



external_function_dims *assign_external_function_dims(int num_inputs, int num_outputs, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    external_function_dims *dims = (external_function_dims *) c_ptr;
    c_ptr += sizeof(external_function_dims);

    dims->input_dims = (int *) c_ptr;
    c_ptr += num_inputs * sizeof(int);

    dims->output_dims = (int *) c_ptr;
    c_ptr += num_outputs * sizeof(int);
    
    assert((char *) raw_memory + external_function_dims_calculate_size(num_inputs, num_outputs) == c_ptr);

    dims->num_inputs = num_inputs;
    dims->num_outputs = num_outputs;

    return dims;
}



int external_function_in_calculate_size(external_function_dims *dims)
{
    int size = sizeof(external_function_in);

    size += dims->num_inputs * sizeof(double *);

    for (int i=0; i<dims->num_inputs; i++) {
        size += dims->input_dims[i] * sizeof(double);
    }

    size += dims->num_outputs * sizeof(bool);

    make_int_multiple_of(8, &size);

    return size;
}



external_function_in *assign_external_function_in(external_function_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    external_function_in *ef_in = (external_function_in *) c_ptr;
    c_ptr += sizeof(external_function_in);

    ef_in->inputs = (double **) c_ptr;
    c_ptr += dims->num_inputs * sizeof(double *);

    for (int i=0; i<dims->num_inputs; i++) {
        ef_in->inputs[i] = (double *) c_ptr;
        c_ptr += dims->input_dims[i] * sizeof(double);
    }

    ef_in->compute_output = (bool *) c_ptr;
    c_ptr += dims->num_outputs * sizeof(bool);
    
    assert((char *) raw_memory + external_function_in_calculate_size(dims) == c_ptr);

    return ef_in;
}



int external_function_out_calculate_size(external_function_dims *dims)
{
    int size = sizeof(external_function_out);

    size += dims->num_outputs * sizeof(double *);

    for (int i=0; i<dims->num_outputs; i++) {
        size += dims->output_dims[i] * sizeof(double);
    }

    make_int_multiple_of(8, &size);

    return size;
}



external_function_out *assign_external_function_out(external_function_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *)raw_memory;

    external_function_out *ef_out = (external_function_out *) c_ptr;
    c_ptr += sizeof(external_function_out);

    ef_out->outputs = (double **) c_ptr;
    c_ptr += dims->num_outputs * sizeof(double *);

    for (int i=0; i<dims->num_outputs; i++) {
        ef_out->outputs[i] = (double *) c_ptr;
        c_ptr += dims->output_dims[i] * sizeof(double);
    }

    assert((char *)raw_memory + external_function_out_calculate_size(dims) == c_ptr);

    return ef_out;
}