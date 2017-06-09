#!/usr/bin/env python

# Tested with both Python 2.7.6 and Python 3.4.3
#
# This Python code collects the source code for testing acados
# on microcontrollers, putting all the necessary C files in
# one directory, and header files in the sub-directory include.
#
# The idea is that when compiling the testing code of acados for
# embedded platforms, when "make" does not fully function like
# on standard Linux platform, all the source code available in
# one directory would allow the compiler to process the code
# easier.
#
# To use for ESP32:
#
# Example usage:
# Assume the source directory of acados is: ~/acados
# The target folder to be created is: mass_spring_qpoases
# This command should be used:
# python test_mass_spring_qpoases.py ~/acados mass_spring_qpoases
#
# Author: Dang Doan
# Date: 2017.04.03-2017.06.09

import sys
import os
import glob

print('Running python script to grab mass_spring_qpoases...')

print(sys.version)  # get python version, for debugging

if len(sys.argv) != 3:
    raise SyntaxError('This script needs exactly 2 arguments: \n \
    test_mass_spring_qpoases.py <acados_top_dir> <new_target_dir>\n \
    Example:\n \
    test_mass_spring_qpoases.py ~/acados mass_spring_qpoases')

# 1. Bring all necessary files to one directory.

top_dir = str(sys.argv[1]).rstrip('/')  # no trailing / in top_dir
target_dir = str(sys.argv[2]).rstrip('/')  # no trailing / in target_dir
# List of file to collect
#  Note: this hard-coded path doesnot work with Windows
workingsourcefiles = [
    'examples/c/mass_spring_qpoases.c',

    'acados/utils/print.c',
    'acados/utils/timing.c',
    'acados/utils/tools.c',
    'acados/ocp_qp/condensing.c',
    'acados/ocp_qp/condensing_helper_functions.c',
    'acados/ocp_qp/ocp_qp_condensing_qpoases.c',

    'external/hpmpc/auxiliary/d_aux_extern_depend_lib4.c',

    'external/blasfeo/auxiliary/i_aux_ext_dep_lib.c',

    'external/qpOASES/src/Constraints.c',
    'external/qpOASES/src/Bounds.c',
    'external/qpOASES/src/Flipper.c',
    'external/qpOASES/src/Indexlist.c',
    'external/qpOASES/src/Matrices.c',
    'external/qpOASES/src/MessageHandling.c',
    'external/qpOASES/src/Options.c',
    'external/qpOASES/src/QProblem.c',
    'external/qpOASES/src/QProblemB.c',
    'external/qpOASES/src/Utils.c'
    ]
workingincludefiles = [
    'acados/ocp_qp/ocp_qp_common.h',
    'acados/ocp_qp/condensing.h',
    'acados/ocp_qp/ocp_qp_condensing_qpoases.h',
    'acados/utils/print.h',
    'acados/utils/types.h',
    'acados/utils/timing.h',
    'acados/utils/tools.h',

    'external/hpmpc/include/aux_d.h',
    'external/hpmpc/include/block_size.h',
    'external/hpmpc/include/kernel_d_lib4.h',

    'external/blasfeo/include/blasfeo_d_aux_ext_dep.h',
    'external/blasfeo/include/blasfeo_i_aux_ext_dep.h',

    'external/qpOASES/include/qpOASES_e/Bounds.h',
    'external/qpOASES/include/qpOASES_e/Constants.h',
    'external/qpOASES/include/qpOASES_e/ConstraintProduct.h',
    'external/qpOASES/include/qpOASES_e/Constraints.h',
    'external/qpOASES/include/qpOASES_e/Flipper.h',
    'external/qpOASES/include/qpOASES_e/Indexlist.h',
    'external/qpOASES/include/qpOASES_e/Matrices.h',
    'external/qpOASES/include/qpOASES_e/MessageHandling.h',
    'external/qpOASES/include/qpOASES_e/Options.h',
    'external/qpOASES/include/qpOASES_e/QProblem.h',
    'external/qpOASES/include/qpOASES_e/QProblemB.h',
    'external/qpOASES/include/qpOASES_e/Utils.h'
    ]
# Files that should be renamed to avoid conflicts
oldfiles = ['external/qpOASES/include/qpOASES_e/Types.h']
newfiles = ['include/qpOASES_e_Types.h']

# Create directory structure and copy files
if not os.path.exists(target_dir):
    os.system('mkdir '+target_dir)
for filename in workingsourcefiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir)
if not os.path.exists(target_dir+'/include'):
    os.system('mkdir '+target_dir+'/include')
for filename in workingincludefiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir+'/include/')
for kk in range(len(oldfiles)):
    os.system('cp '+top_dir+'/'+oldfiles[kk]+' '+target_dir+'/'+newfiles[kk])

print('Step 1: Necessary files copied.')

# 2. Modify .h and .c files to adapt to the new code structure:
# List of texts to be replaced:
old_text = [
    'acados/ocp_qp/condensing.h',
    'acados/ocp_qp/condensing_helper_functions.c',
    'acados/ocp_qp/ocp_qp_common.h',
    'acados/ocp_qp/ocp_qp_condensing_qpoases.h',
    'acados/ocp_qp/ocp_qp_hpmpc.h',
    'acados/utils/print.h',
    'acados/utils/timing.h',
    'acados/utils/tools.h',
    'acados/utils/types.h',

    'hpmpc/include/aux_d.h',
    '../include/block_size.h',
    '../include/kernel_d_lib4.h',

    'blasfeo/include/blasfeo_common.h',
    'blasfeo/include/blasfeo_d_aux_ext_dep.h',
    'blasfeo/include/blasfeo_i_aux_ext_dep.h',

    'qpOASES_e/Bounds.h',
    'qpOASES_e/Constants.h',
    'qpOASES_e/Constraints.h',
    'qpOASES_e/ConstraintProduct.h',
    'qpOASES_e/Flipper.h',
    'qpOASES_e/Indexlist.h',
    'qpOASES_e/Matrices.h',
    'qpOASES_e/MessageHandling.h',
    'qpOASES_e/Options.h',
    'qpOASES_e/QProblem.h',
    'qpOASES_e/QProblemB.h',
    'qpOASES_e/Types.h',
    'qpOASES_e/Utils.h'
    ]
# List of new texts to replace old ones,
#  in corresponding order to old_text:
new_text = [
    'condensing.h',
    'condensing_helper_functions.c',
    'ocp_qp_common.h',
    'ocp_qp_condensing_qpoases.h',
    'ocp_qp_hpmpc.h',
    'print.h',
    'timing.h',
    'tools.h',
    'types.h',

    'aux_d.h',
    'block_size.h',
    'kernel_d_lib4.h',

    'blasfeo_common.h',
    'blasfeo_d_aux_ext_dep.h',
    'blasfeo_i_aux_ext_dep.h',

    'Bounds.h',
    'Constants.h',
    'Constraints.h',
    'ConstraintProduct.h',
    'Flipper.h',
    'Indexlist.h',
    'Matrices.h',
    'MessageHandling.h',
    'Options.h',
    'QProblem.h',
    'QProblemB.h',
    'qpOASES_e_Types.h',  # This filename is changed
    'Utils.h'
    ]

len_old_text = len(old_text)
len_new_text = len(new_text)

if len_old_text != len_new_text:
    raise ValueError('Number of old and new texts not match')

files = glob.glob(target_dir+"/*.c")
for file in files:
    with open(file) as objFile:
        txtFile = objFile.read()
    for replacetext in range(len_old_text):
        txtFile = txtFile.replace(old_text[replacetext], new_text[replacetext])
    with open(file, "w") as objFile:
        objFile.write(txtFile)

files = glob.glob(target_dir+"/include/*.h")
for file in files:
    with open(file) as objFile:
        txtFile = objFile.read()
    for replacetext in range(len_old_text):
        txtFile = txtFile.replace(old_text[replacetext], new_text[replacetext])
    with open(file, "w") as objFile:
        objFile.write(txtFile)

print('Step 2: Path information in files modified to the new structure.')

# 3. Add specific code to HPMPC and BLASFEO files:
# List of files to be modified:
files = [
    'include/block_size.h',
    'include/Constants.h',
    'include/ocp_qp_condensing_qpoases.h'
    ]
# List of lines to be added in the beginning of files,
#  in corresponding order with the list files:
lines = [
    '#include "target.h"\n',
    '#include "ocp_qp_condensing_qpoases.h"\n',
    '#define QPOASES_NVMAX 60\n#define QPOASES_NCMAX 228\n'
    ]

if len(files) != len(lines):
    raise ValueError('Number of files and added lines not match')

for kk in range(len(files)):
    with open(target_dir+'/'+files[kk]) as objFile:
        txtFile = objFile.read()
    with open(target_dir+'/'+files[kk], "w") as objFile:
        objFile.write(lines[kk])  # write the line to the beginning
        objFile.write(txtFile)

print('Step 3: Common header file included in specific files.')

# 4. Copy Makefile and specific setting files
os.system('cp '+top_dir+'/experimental/dang/esp32/script_mass_spring_qpoases/Makefile '+target_dir)
os.system('cp '+top_dir+'/experimental/dang/esp32/script_mass_spring_qpoases/target.h '+target_dir+'/include/')

print('Step 4: Makefile, and HPMPC target.h replaced.')

# 5. Display further instructions
print('Please do next steps in terminal:')
print(' cd '+target_dir)
print(' make')
print('Then run the binary file in '+target_dir+'/bin')
print('To remove binary objects: make clean\n')
