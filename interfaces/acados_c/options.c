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

#include "acados_c/options.h"

#include <string.h>



int get_option_int(const void *args_, const char *option)
{
    return 0;
}



int set_option_int(void *args_, const char *option, const int value)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}



const int *get_option_int_array(const void *args_, const char *option)
{
    return NULL;
}



int set_option_int_array(void *args_, const char *option, const int *value)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}



double get_option_double(const void *args_, const char *option)
{
    return 0;
}



int set_option_double(void *args_, const char *option, const double value)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}



const double *get_option_double_array(const void *args_, const char *option)
{
    return NULL;
}



int set_option_double_array(void *args_, const char *option, const double *value)
{
    int return_value = ACADOS_SUCCESS;

    return return_value;
}
