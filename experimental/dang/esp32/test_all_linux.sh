#!/usr/bin/env bash

PROJECT=chen_nmpc_qpoases

python3 script_chen_nmpc_qpoases/test_chen_nmpc_qpoases.py ../../.. $PROJECT
cd $PROJECT
make "PROJECT=$PROJECT"
bin/$PROJECT.exe
make clean "PROJECT=$PROJECT"
cd ..

PROJECT=mass_spring_hpmpc

python3 script_mass_spring_hpmpc/test_mass_spring_hpmpc.py ../../.. $PROJECT
cd $PROJECT
make "PROJECT=$PROJECT"
bin/$PROJECT.exe
make clean "PROJECT=$PROJECT"
cd ..

PROJECT=mass_spring_qpoases

python3 script_mass_spring_qpoases/test_mass_spring_qpoases.py ../../.. $PROJECT
cd $PROJECT
make "PROJECT=$PROJECT"
bin/$PROJECT.exe
make clean "PROJECT=$PROJECT"
