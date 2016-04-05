rm *.o *_wrap.c *.so

swig -python erk_integrator.i
gcc -c -fPIC -O3 model.c auxiliary_functions.c timing_functions.c
gcc -c -fPIC -O3 erk_integrator.c erk_integrator_wrap.c -I. -I/usr/local/include/python2.7
gcc -v -undefined dynamic_lookup -dynamiclib model.o erk_integrator.o erk_integrator_wrap.o auxiliary_functions.o timing_functions.o -o _erk_integrator.so
