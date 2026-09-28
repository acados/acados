/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * Licensed under the 2-Clause BSD License.
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
