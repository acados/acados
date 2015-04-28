rm *.o *_wrap.c *.so

swig -python sparse_erk_integrator.i
swig -python erk_integrator.i
swig -python rk4_integrator.i
gcc -c -fPIC -O3 model.c auxiliary_functions.c timing_functions.c
gcc -c -fPIC -O3 erk_integrator.c erk_integrator_wrap.c -I. -I/usr/include/python2.7
gcc -c -fPIC -O3 sparse_erk_integrator.c sparse_erk_integrator_wrap.c -I. -I/usr/include/python2.7
gcc -c -fPIC -O3 rk4_integrator.c rk4_integrator_wrap.c -I. -I/usr/include/python2.7
ld -shared model.o erk_integrator.o erk_integrator_wrap.o auxiliary_functions.o timing_functions.o -o _erk_integrator.so 
ld -shared model.o sparse_erk_integrator.o sparse_erk_integrator_wrap.o auxiliary_functions.o timing_functions.o -o _sparse_erk_integrator.so 
ld -shared model.o rk4_integrator.o rk4_integrator_wrap.o auxiliary_functions.o timing_functions.o -o _rk4_integrator.so 

