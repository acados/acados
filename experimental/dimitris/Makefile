ACADOS=$(abspath ../../)
ACADOS_EXT=$(ACADOS)/external

OOQP=$(ACADOS_EXT)/OOQP
OOQPINCLUDEDIR=$(OOQP)/include
OOQPLIBDIR=$(OOQP)/lib
SORTOBJ = $(OOQP)/src/Utilities/sort.o

BLASFEO=$(ACADOS_EXT)/blasfeo
BLASFEOLIBDIR=$(BLASFEO)

HPMPC=$(ACADOS_EXT)/hpmpc
HPMPCLIBDIR=$(HPMPC)

QPOASES=$(ACADOS_EXT)/qpOASES
QPOASESINCLUDEDIR=$(QPOASES)/include
QPOASESLIBDIR=$(QPOASES)/bin

LDFLAGS =-L$(OOQPLIBDIR)
LDFLAGS+=-L$(BLASFEOLIBDIR)
LDFLAGS+=-L$(QPOASESLIBDIR)
LDFLAGS+=-L$(HPMPCLIBDIR)

BLAS     = -lblas
MA27LIB  = $(OOQP)/libma27.a
FLIBS    = -lgfortran
# FLIBS = /usr/lib/gcc/x86_64-linux-gnu/5/libgfortran.so
OOQP_SPARSE_LIBS = -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase \
                 $(BLAS) $(MA27LIB) $(FLIBS)

LIBS= -lqpOASES_e -lhpmpc -lblasfeo
LIBS+=$(OOQP_SPARSE_LIBS)
LIBS+=-lm

toy_memory_allocation: toy_memory_allocation.c
	gcc -g -O0 toy_memory_allocation.c -o toy_memory_allocation.exe
	valgrind --leak-check=full ./toy_memory_allocation.exe

ooqp: test_qp.c
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_ooqp.c -o ocp_qp_ooqp.o -I$(OOQPINCLUDEDIR) -I$(ACADOS)  -I$(ACADOS_EXT)
	g++ -c -g ../../acados/ocp_qp/ocp_qp_ooqp_helper.cpp -o ocp_qp_ooqp_helper.o -I$(ACADOS_EXT) -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=1 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	g++ -g -O0 ocp_qp_ooqp.o ocp_qp_ooqp_helper.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o $(SORTOBJ) -o test_qp.exe $(LDFLAGS) $(LIBS)
	valgrind --leak-check=full ./test_qp.exe

qpoases: test_qp.c
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_condensing_qpoases.c -o ocp_qp_condensing_qpoases.o -I$(ACADOS_EXT) -I$(ACADOS) -I$(QPOASESINCLUDEDIR)
	gcc -c -g -O0 ../../acados/ocp_qp/condensing.c -o ocp_qp_condensing.o -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=2 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -g -O0 ocp_qp_condensing.o ocp_qp_condensing_qpoases.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o -o test_qp.exe $(LDFLAGS) $(LIBS)
	valgrind --leak-check=full ./test_qp.exe

hpmpc: test_qp.c # TODO FIX!
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_hpmpc.c -o ocp_qp_hpmpc.o -I$(ACADOS_EXT) -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=3 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -g -O0 ocp_qp_hpmpc.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o -o test_qp.exe $(LDFLAGS) $(LIBS)
	./test_qp.exe

clean:
	rm *.exe *.o
