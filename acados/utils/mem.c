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

#include <stdio.h>

#include "acados/utils/mem.h"


void make_int_multiple_of(int num, int *size) {
    *size = (*size + num - 1) / num * num;
}



void align_char_to(int num, char **c_ptr) {
    size_t s_ptr = (size_t)*c_ptr;
    s_ptr = (s_ptr + num - 1) / num * num;
    *c_ptr = (char *)s_ptr;
}