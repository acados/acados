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
# The target folder to be created is: ocp_qp_hpmpc
# This command should be used:
# python test_ocp_qp_hpmpc.py ~/acados ocp_qp_hpmpc
#
# Author: Dang Doan
# Date: 2017.04.03

import sys
import os
import glob
from subprocess import call
from os.path import join

print('Running python script to grab ocp_qp_hpmpc...')

print(sys.version) # get python version, for debugging

if len(sys.argv)!= 3:
    raise SyntaxError('This script needs exactly 2 arguments: \n \
    test_ocp_qp_hpmpc.py <acados_top_dir> <new_target_dir>\n \
    Example:\n \
    test_ocp_qp_hpmpc.py ~/acados ocp_qp_hpmpc')

# 1. Bring all necessary files to one directory.

top_dir = str(sys.argv[1]).rstrip('/') # no trailing / in top_dir
target_dir = str(sys.argv[2]).rstrip('/') # no trailing / in target_dir
# List of file to collect
#  Note: this hard-coded path doesnot work with Windows
workingcodefiles = \
    ['examples/c/mass_spring_hpmpc.c', \
    'acados/ocp_qp/ocp_qp_hpmpc.c', \
    'acados/utils/print.c', \
    'acados/utils/timing.c', \
    'acados/utils/math.c', \
    'external/hpmpc/interfaces/c/fortran_order_interface_libstr.c', \
    'external/hpmpc/mpc_solvers/d_ip2_res_hard_libstr.c', \
    'external/hpmpc/mpc_solvers/c99/d_res_ip_res_hard_libstr.c', \
    'external/hpmpc/mpc_solvers/c99/d_aux_ip_hard_libstr.c', \
    'external/hpmpc/lqcp_solvers/d_back_ric_rec_libstr.c', \
    'external/hpmpc/lqcp_solvers/d_part_cond_libstr.c', \
    'external/hpmpc/auxiliary/d_aux_extern_depend_lib4.c', \
    'external/hpmpc/auxiliary/i_aux.c', \
    'external/blasfeo/kernel/c99/kernel_sgemm_4x4_lib4.c', \
    'external/blasfeo/kernel/c99/kernel_sgemv_4_lib4.c', \
    'external/blasfeo/blas/d_blas1_lib.c', \
    'external/blasfeo/blas/d_blas2_lib.c', \
    'external/blasfeo/blas/d_blas3_diag_lib.c', \
    'external/blasfeo/blas/d_blas3_lib.c', \
    'external/blasfeo/blas/d_lapack_lib.c', \
    'external/blasfeo/auxiliary/d_aux_lib.c', \
    'external/blasfeo/auxiliary/d_aux_extern_depend_lib.c']
workingheaderfiles =\
    ['acados/ocp_qp/ocp_qp_common.h', \
    'acados/ocp_qp/ocp_qp_hpmpc.h', \
    'acados/sim/sim_rk_common.h', \
    'acados/utils/types.h', \
    'acados/utils/print.h', \
    'acados/utils/timing.h', \
    'acados/utils/math.h', \
    'external/hpmpc/include/aux_s.h', \
    'external/hpmpc/include/blas_d.h', \
    'external/hpmpc/include/block_size.h', \
    'external/hpmpc/include/target.h', \
    'external/hpmpc/include/c_interface.h', \
    'external/hpmpc/include/d_blas_aux.h', \
    'external/hpmpc/include/tree.h', \
    'external/hpmpc/include/kernel_d_lib4.h', \
    'external/hpmpc/include/lqcp_aux.h', \
    'external/hpmpc/include/lqcp_solvers.h', \
    'external/hpmpc/include/mpc_aux.h', \
    'external/hpmpc/include/mpc_solvers.h', \
    'external/hpmpc/include/aux_d.h', \
    'external/blasfeo/include/blasfeo_common.h', \
    'external/blasfeo/include/blasfeo_target.h', \
    'external/blasfeo/include/blasfeo_block_size.h', \
    'external/blasfeo/include/blasfeo_d_kernel.h', \
    'external/blasfeo/include/blasfeo_d_blas.h', \
    'external/blasfeo/include/blasfeo_d_aux.h', \
    'external/blasfeo/include/blasfeo_i_aux.h']
# Create directory structure and copy files
if not os.path.exists(target_dir):
    os.system('mkdir '+target_dir)
for filename in workingcodefiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir)
if not os.path.exists(target_dir+'/include'):
    os.system('mkdir '+target_dir+'/include')
for filename in workingheaderfiles:
    os.system('cp '+top_dir+'/'+filename+' '+target_dir+'/include/')

print('Step 1: Necessary files copied.')

# 2. Modify .h and .c files to adapt to the new code structure:
# List of texts to be replaced:
old_text =\
    ['acados/ocp_qp/ocp_qp_common.h', \
    'acados/ocp_qp/ocp_qp_hpmpc.h', \
    'acados/utils/print.h', \
    'acados/utils/timing.h', \
    'acados/utils/math.h', \
    'acados/utils/types.h', \
    'acados/sim/sim_common.h', \

    'hpmpc/include/aux_d.h', \
    'hpmpc/include/c_interface.h', \
    'hpmpc/include/lqcp_solvers.h', \
    'hpmpc/include/mpc_aux.h', \
    'hpmpc/include/mpc_solvers.h', \
    '../../include/aux_d.h', \
    '../../include/aux_s.h', \
    '../../include/blas_d.h', \
    '../../include/block_size.h', \
    '../../include/lqcp_solvers.h', \
    '../../include/mpc_aux.h', \
    '../../include/mpc_solvers.h', \
    '../../include/target.h', \
    '../include/aux_d.h', \
    '../include/blas_d.h', \
    '../include/block_size.h', \
    '../include/d_blas_aux.h', \
    '../include/kernel_d_lib4.h', \
    '../include/lqcp_aux.h', \
    '../include/lqcp_solvers.h', \
    '../include/mpc_aux.h', \

    'blasfeo/include/blasfeo_common.h', \
    'blasfeo/include/blasfeo_d_aux.h', \
    'blasfeo/include/blasfeo_d_blas.h', \
    'blasfeo/include/blasfeo_target.h', \
    '../include/blasfeo_block_size.h', \
    '../include/blasfeo_common.h', \
    '../include/blasfeo_d_aux.h', \
    '../include/blasfeo_d_kernel.h' \
    ]
# List of new texts to replace old ones,
#  in corresponding order to old_text:
new_text =\
    ['ocp_qp_common.h', \
    'ocp_qp_hpmpc.h', \
    'print.h', \
    'timing.h', \
    'tools.h', \
    'types.h', \
    'sim_common.h', \

    'aux_d.h', \
    'c_interface.h', \
    'lqcp_solvers.h', \
    'mpc_aux.h', \
    'mpc_solvers.h', \
    'aux_d.h', \
    'aux_s.h', \
    'blas_d.h', \
    'block_size.h', \
    'lqcp_solvers.h', \
    'mpc_aux.h', \
    'mpc_solvers.h', \
    'target.h', \
    'aux_d.h', \
    'blas_d.h', \
    'block_size.h', \
    'd_blas_aux.h', \
    'kernel_d_lib4.h', \
    'lqcp_aux.h', \
    'lqcp_solvers.h', \
    'mpc_aux.h', \

    'blasfeo_common.h', \
    'blasfeo_d_aux.h', \
    'blasfeo_d_blas.h', \
    'blasfeo_target.h', \
    'blasfeo_block_size.h', \
    'blasfeo_common.h', \
    'blasfeo_d_aux.h', \
    'blasfeo_d_kernel.h' \
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
files =\
    ['include/blasfeo_common.h', \
    'include/blas_d.h', \
    'include/aux_d.h', \
    'include/block_size.h', \
    'fortran_order_interface_libstr.c', \
    'd_ip2_res_hard_libstr.c', \
    'd_res_ip_res_hard_libstr.c', \
    'd_back_ric_rec_libstr.c', \
    'd_aux_ip_hard_libstr.c']
# List of lines to be added in the beginning of files,
#  in corresponding order with the list files:
lines =\
    ['#include "blasfeo_target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n', \
    '#include "target.h"\n']

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
os.system('cp '+top_dir+'/experimental/dang/esp32/script_test_ocp_qp_hpmpc/Makefile '+target_dir)
os.system('cp '+top_dir+'/experimental/dang/esp32/script_test_ocp_qp_hpmpc/target.h '+target_dir+'/include/')
os.system('cp '+top_dir+'/experimental/dang/esp32/script_test_ocp_qp_hpmpc/blasfeo_target.h '+target_dir+'/include/')

print('Step 4: Makefile, blasfeo_target.h and HPMPC target.h replaced.')

# 5. Display further instruction
print('Please do next steps in terminal:')
print(' cd '+target_dir)
print(' make')
print('Then run the binary file in '+target_dir+'/bin')
print('To remove binary objects: make clean\n')
