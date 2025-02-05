/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#ifndef ACADOS_UTILS_STRSEP_H_
#define ACADOS_UTILS_STRSEP_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>

// Inline function definition
static inline void extract_module_name(const char *module_field, char *module, int *module_length, char **ptr_module)
{
    // extract module name from string of the form <module>_<field>
    char *char_ = strchr(module_field, '_');
    if (char_ != NULL)
    {
        *module_length = char_ - module_field;
        // Copy the module name into the module array
        strncpy(module, module_field, *module_length);
        module[*module_length] = '\0'; // add end of string
        *ptr_module = module;
    }
}


static inline void extract_field_name(const char *module_field, char *field, int *field_length, char **ptr_field)
{
    // extract field name from string of the form <module>_<field>
    char *char_ = strchr(module_field, '_');
    if (char_ != NULL)
    {
        int length_prefix = char_ - module_field + 1;
        *field_length = strlen(module_field) - length_prefix;
        // Copy the field name into the module array
        strncpy(field, module_field+length_prefix, *field_length);
        field[*field_length] = '\0'; // add end of string
        *ptr_field = field;
    }
}

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_UTILS_STRSEP_H_
