#!/bin/bash
python3 script_test_nmpc/test_nmpc.py ../../.. test_nmpc
cd test_nmpc
make
bin/test_nmpc.exe
make clean
cd ..

python3 script_test_ocp_qp_hpmpc/test_ocp_qp_hpmpc.py ../../.. test_ocp_qp_hpmpc
cd test_ocp_qp_hpmpc
make
bin/test_ocp_qp_hpmpc.exe
make clean
