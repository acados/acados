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
# The target folder to be created is: chen_nmpc_qpoases
# This command should be used:
# python test_nmpc_qpoases.py ~/acados chen_nmpc_qpoases
#
# Author: Dang Doan
# Date: 2017.04.03

import sys
import os
import glob
from subprocess import call
from os.path import join

print('Running python script to grab chen_nmpc_qpoases...')

print(sys.version) # get python version, for debugging

if len(sys.argv)!= 3:
    raise SyntaxError('This script needs exactly 2 arguments: \n \
    test_nmpc_qpoases.py <acados_top_dir> <new_target_dir>\n \
    Example:\n \
    test_nmpc_qpoases.py ~/acados chen_nmpc_qpoases')

# 1. Bring all necessary files to one directory.

top_dir = str(sys.argv[1]).rstrip('/') # no trailing / in top_dir
target_dir = str(sys.argv[2]).rstrip('/') # no trailing / in target_dir
# List of file to collect
#  Note: this hard-coded path doesnot work with Windows
workingcodefiles = [\
    'examples/c/chen_nmpc_qpoases.c', \
    'examples/c/Chen_model/chen_model.c', \

    'acados/utils/print.c', \
    'acados/utils/timing.c', \
    'acados/ocp_qp/condensing.c', \
    'acados/ocp_qp/condensing_helper_functions.c', \
    'acados/ocp_qp/ocp_qp_condensing_qpoases.c', \
    'acados/sim/sim_erk_integrator.c', \

    'external/hpmpc/auxiliary/d_aux_extern_depend_lib4.c', \

    'external/blasfeo/auxiliary/i_aux_extern_depend_lib.c', \

    'external/qpOASES/src/Constraints.c', \
    'external/qpOASES/src/Bounds.c', \
    'external/qpOASES/src/Flipper.c', \
    'external/qpOASES/src/Indexlist.c', \
    'external/qpOASES/src/Matrices.c', \
    'external/qpOASES/src/MessageHandling.c', \
    'external/qpOASES/src/Options.c', \
    'external/qpOASES/src/QProblem.c', \
    'external/qpOASES/src/QProblemB.c', \
    'external/qpOASES/src/Utils.c' \
    ]
workingheaderfiles = [\
    'examples/c/Chen_model/chen_model.h', \
    'acados/ocp_qp/ocp_qp_common.h', \
    'acados/ocp_qp/condensing.h', \
    'acados/ocp_qp/ocp_qp_condensing_qpoases.h', \
    'acados/sim/sim_common.h', \
    'acados/sim/sim_erk_integrator.h', \
    'acados/sim/sim_collocation.h', \
    'acados/sim/sim_rk_common.h', \
    'acados/utils/print.h', \
    'acados/utils/types.h', \
    'acados/utils/timing.h', \

    'external/hpmpc/include/aux_d.h', \
    'external/hpmpc/include/block_size.h', \
    'external/hpmpc/include/kernel_d_lib4.h', \

    'external/blasfeo/include/blasfeo_i_aux.h', \

    'external/qpOASES/include/qpOASES_e/Bounds.h', \
    'external/qpOASES/include/qpOASES_e/Constants.h', \
    'external/qpOASES/include/qpOASES_e/ConstraintProduct.h', \
    'external/qpOASES/include/qpOASES_e/Constraints.h', \
    'external/qpOASES/include/qpOASES_e/Flipper.h', \
    'external/qpOASES/include/qpOASES_e/Indexlist.h', \
    'external/qpOASES/include/qpOASES_e/Matrices.h', \
    'external/qpOASES/include/qpOASES_e/MessageHandling.h', \
    'external/qpOASES/include/qpOASES_e/Options.h', \
    'external/qpOASES/include/qpOASES_e/QProblem.h', \
    'external/qpOASES/include/qpOASES_e/QProblemB.h', \
    'external/qpOASES/include/qpOASES_e/Utils.h' \
    ]
# Files that should be renamed to avoid conflicts
oldfiles = ['external/qpOASES/include/qpOASES_e/Types.h']
newfiles = ['include/qpOASES_e_Types.h']

# Create directory structure and copy files
if not os.path.exists(target_dir):
    os.system('mkdir '+target_dir)
for filename in workingcodefiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir)
if not os.path.exists(target_dir+'/include'):
    os.system('mkdir '+target_dir+'/include')
for filename in workingheaderfiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir+'/include/')
for kk in range(len(oldfiles)):
    os.system('cp '+top_dir+'/'+oldfiles[kk]+' '+target_dir+'/'+newfiles[kk])

print('Step 1: Necessary files copied.')

# 2. Modify .h and .c files to adapt to the new code structure:
# List of texts to be replaced:
old_text = [\
    'examples/c/Chen_model/chen_model.h', \
    'acados/ocp_qp/condensing.h', \
    'acados/ocp_qp/condensing_helper_functions.c', \
    'acados/ocp_qp/ocp_qp_common.h', \
    'acados/ocp_qp/ocp_qp_condensing_qpoases.h', \
    'acados/ocp_qp/ocp_qp_hpmpc.h', \
    'acados/sim/sim_common.h', \
    'acados/sim/sim_erk_integrator.h', \
    'acados/sim/sim_collocation.h', \
    'acados/sim/sim_rk_common.h', \
    'acados/utils/print.h', \
    'acados/utils/timing.h', \
    'acados/utils/types.h', \

    'hpmpc/include/aux_d.h', \
    '../include/block_size.h', \
    '../include/kernel_d_lib4.h', \

    'blasfeo/include/blasfeo_common.h', \
    'blasfeo/include/blasfeo_i_aux.h', \

    'qpOASES_e/Bounds.h', \
    'qpOASES_e/Constants.h', \
    'qpOASES_e/Constraints.h', \
    'qpOASES_e/ConstraintProduct.h', \
    'qpOASES_e/Flipper.h', \
    'qpOASES_e/Indexlist.h', \
    'qpOASES_e/Matrices.h', \
    'qpOASES_e/MessageHandling.h', \
    'qpOASES_e/Options.h', \
    'qpOASES_e/QProblem.h', \
    'qpOASES_e/QProblemB.h', \
    'qpOASES_e/Types.h', \
    'qpOASES_e/Utils.h' \
    ]
# List of new texts to replace old ones,
#  in corresponding order to old_text:
new_text = [\
    'chen_model.h', \
    'condensing.h', \
    'condensing_helper_functions.c', \
    'ocp_qp_common.h', \
    'ocp_qp_condensing_qpoases.h', \
    'ocp_qp_hpmpc.h', \
    'sim_common.h', \
    'sim_erk_integrator.h', \
    'sim_collocation.h', \
    'sim_rk_common.h', \
    'print.h', \
    'timing.h', \
    'types.h', \

    'aux_d.h', \
    'block_size.h', \
    'kernel_d_lib4.h', \

    'blasfeo_common.h', \
    'blasfeo_i_aux.h', \

    'Bounds.h', \
    'Constants.h', \
    'Constraints.h', \
    'ConstraintProduct.h', \
    'Flipper.h', \
    'Indexlist.h', \
    'Matrices.h', \
    'MessageHandling.h', \
    'Options.h', \
    'QProblem.h', \
    'QProblemB.h', \
    'qpOASES_e_Types.h', \
    'Utils.h' \
    ]

len_old_text = len(old_text)
len_new_text = len(new_text)

if len_old_text != len_new_text:
    raise ValueError('Number of old and new texts not match')

files = glob.glob(target_dir+"/*.c")
for file in files:
    objFile = open(file, "r")
    txtFile = objFile.read()
    objFile.close()
    for replacetext in range(len_old_text):
        txtFile = txtFile.replace(old_text[replacetext],new_text[replacetext])
    objFile = open(file, "w")
    objFile.write(txtFile)
    objFile.close()

files = glob.glob(target_dir+"/include/*.h")
for file in files:
    objFile = open(file, "r")
    txtFile = objFile.read()
    objFile.close()
    for replacetext in range(len_old_text):
        txtFile = txtFile.replace(old_text[replacetext],new_text[replacetext])
    objFile = open(file, "w")
    objFile.write(txtFile)
    objFile.close()

print('Step 2: Path information in files modified to the new structure.')

# 3. Add specific code to HPMPC and BLASFEO files:
# List of files to be modified:
files = ['include/block_size.h']
# List of lines to be added in the beginning of files,
#  in corresponding order with the list files:
lines = ['#include "target.h"\n']

if len(files) != len(lines):
    raise ValueError('Number of files and added lines not match')

for kk in range(len(files)):
    objFile = open(target_dir+'/'+files[kk], "r")
    txtFile = objFile.read()
    objFile.close()
    objFile = open(target_dir+'/'+files[kk], "w")
    objFile.write(lines[kk]) # write the line to the beginning
    objFile.write(txtFile)
    objFile.close()

print('Step 3: Common header file included in specific files.')

# 4. Copy Makefile and specific setting files
os.system('cp '+top_dir+'/experimental/dang/esp32/script_test_nmpc_qpoases/Makefile '+target_dir)
os.system('cp '+top_dir+'/experimental/dang/esp32/script_test_nmpc_qpoases/target.h '+target_dir+'/include/')

print('Step 4: Makefile, and HPMPC target.h replaced.')

# 5. Display further instructions
print('Please do next steps in terminal:')
print(' cd '+target_dir)
print(' make')
print('Then run the binary file in '+target_dir+'/bin')
print('To remove binary objects: make clean\n')
