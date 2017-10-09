OS := $(shell uname)

ACADOS=$(abspath ../../)
ACADOS_EXT=$(ACADOS)/external
ACADOSLIBDIR=$(ACADOS)/build/acados

OOQP=$(ACADOS_EXT)/OOQP
OOQPINCLUDEDIR=$(OOQP)/include
OOQPLIBDIR=$(OOQP)/lib
SORTOBJ = $(OOQP)/src/Utilities/sort.o

BLASFEO=$(ACADOS_EXT)/blasfeo
BLASFEOLIBDIR=$(BLASFEO)

HPMPC=$(ACADOS_EXT)/hpmpc
HPMPCLIBDIR=$(HPMPC)

EIGEN=$(ACADOS_EXT)/eigen

QPOASES=$(ACADOS_EXT)/qpOASES
QPOASESINCLUDEDIR=$(QPOASES)/include
QPOASESLIBDIR=$(QPOASES)/bin

QPDUNES=$(ACADOS_EXT)/qpDUNES-dev
QPDUNESINCLUDEDIR=$(QPDUNES)/include
QPDUNESLIBDIR=$(QPDUNES)/build/lib

LDFLAGS =-L$(OOQPLIBDIR)
LDFLAGS+=-L$(BLASFEOLIBDIR)
LDFLAGS+=-L$(QPOASESLIBDIR)
LDFLAGS+=-L$(HPMPCLIBDIR)
LDFLAGS+=-L$(QPDUNESLIBDIR)
LDFLAGS+=-L$(ACADOSLIBDIR)

BLAS     = -lblas
MA27LIB  = $(OOQP)/libma27.a
ifeq ($(OS), Darwin)
	FLIBS = /usr/local/lib/gcc/7/libgfortran.dylib
else
	FLIBS = -lgfortran
endif

OOQP_SPARSE_LIBS = -looqpgensparse -looqpsparse  -looqpgondzio -looqpbase \
                 $(BLAS) $(MA27LIB) $(FLIBS)

LIBS= -lqpOASES_e -lhpmpc -lblasfeo -lqpdunes -lacados
LIBS+=$(OOQP_SPARSE_LIBS)
LIBS+=-lm

toy_memory_allocation: toy_memory_allocation.c
	gcc -g -O0 toy_memory_allocation.c -o toy_memory_allocation.exe
	valgrind --leak-check=full ./toy_memory_allocation.exe

qpdunes: test_qpdunes.c
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_qpdunes.c -o ocp_qp_qpdunes.o -I$(ACADOS) -I$(ACADOS_EXT) -I$(QPDUNESINCLUDEDIR)
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_ooqp.c -o ocp_qp_ooqp.o -I$(OOQPINCLUDEDIR) -I$(ACADOS)  -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../acados/utils/timing.c -o timing.o -I$(ACADOS)
	gcc -c -g -O0 test_qpdunes.c -o test_qpdunes.o -I$(ACADOS) -I$(ACADOS_EXT) -I$(QPDUNESINCLUDEDIR)
	g++ -g -O0 timing.o ocp_qp_qpdunes.o ocp_qp_ooqp.o allocate_ocp_qp.o read_ocp_qp_in.o test_qpdunes.o -o qpdunes.exe $(LDFLAGS) $(LIBS)
	./qpdunes.exe
	# valgrind --leak-check=full ./qpdunes.exe

ooqp: test_qp.c
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_ooqp.c -o ocp_qp_ooqp.o -I$(OOQPINCLUDEDIR) -I$(ACADOS)  -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../acados/utils/timing.c -o timing.o -I$(ACADOS)
	gcc -c -g -O0 ../../acados/ocp_qp/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=1 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	g++ -g -O0 ocp_qp_ooqp.o timing.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o $(SORTOBJ) -o test_qp.exe $(LDFLAGS) $(LIBS)
	# valgrind --leak-check=full ./test_qp.exe
	./test_qp.exe

qpoases: test_qp.c
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_condensing_qpoases.c -o ocp_qp_condensing_qpoases.o -I$(ACADOS_EXT) -I$(ACADOS) -I$(QPOASESINCLUDEDIR)
	gcc -c -g -O0 ../../acados/ocp_qp/condensing.c -o ocp_qp_condensing.o -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=2 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -g -O0 ocp_qp_condensing.o ocp_qp_condensing_qpoases.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o -o test_qp.exe $(LDFLAGS) $(LIBS)
	valgrind --leak-check=full ./test_qp.exe

hpmpc: test_qp.c # TODO FIX FOR MAC!
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_hpmpc.c -o ocp_qp_hpmpc.o -I$(ACADOS_EXT) -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 -DSOLVER=3 test_qp.c -o test_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -g -O0 ocp_qp_hpmpc.o allocate_ocp_qp.o read_ocp_qp_in.o test_qp.o -o test_qp.exe $(LDFLAGS) $(LIBS)
	./test_qp.exe

all: test_all_qp_solvers.c
	gcc -c -g -O0 ../../acados/utils/timing.c -o timing.o -I$(ACADOS)
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_ooqp.c -o ocp_qp_ooqp.o -I$(OOQPINCLUDEDIR) -I$(ACADOS)  -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_condensing_qpoases.c -o ocp_qp_condensing_qpoases.o -I$(ACADOS_EXT) -I$(ACADOS) -I$(QPOASESINCLUDEDIR)
	gcc -c -g -O0 ../../acados/ocp_qp/condensing.c -o ocp_qp_condensing.o -I$(ACADOS)
	gcc -c -g -O0 ../../acados/ocp_qp/ocp_qp_hpmpc.c -o ocp_qp_hpmpc.o -I$(ACADOS_EXT) -I$(ACADOS)
	gcc -c -g -O0 ../../acados/utils/allocate_ocp_qp.c -o allocate_ocp_qp.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 test_all_qp_solvers.c -o test_all_qp_solvers.o -I$(ACADOS) -I$(ACADOS_EXT)
	g++ -g -O0 timing.o ocp_qp_ooqp.o ocp_qp_condensing.o ocp_qp_condensing_qpoases.o ocp_qp_hpmpc.o allocate_ocp_qp.o read_ocp_qp_in.o test_all_qp_solvers.o $(SORTOBJ) -o test_qp.exe $(LDFLAGS) $(LIBS)
	# valgrind --leak-check=full ./test_qp.exe
	./test_qp.exe

chain:
	gcc -c -g -O0 $(ACADOS)/test/ocp_nlp/chain/Chain_model.c -o model.o -I$(ACADOS)
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	gcc -c -g -O0 $(ACADOS)/build/test/jac_chain_nm2.c -o jac_chain_nm2.o
	gcc -c -g -O0 $(ACADOS)/build/test/jac_chain_nm3.c -o jac_chain_nm3.o
	gcc -c -g -O0 $(ACADOS)/build/test/jac_chain_nm4.c -o jac_chain_nm4.o
	gcc -c -g -O0 $(ACADOS)/build/test/vde_chain_nm2.c -o vde_chain_nm2.o
	gcc -c -g -O0 $(ACADOS)/build/test/vde_chain_nm3.c -o vde_chain_nm3.o
	gcc -c -g -O0 $(ACADOS)/build/test/vde_chain_nm4.c -o vde_chain_nm4.o
	g++ -c -g -O0 read_matrix_mod.cpp -o read_matrix.o -I$(EIGEN) -I$(ACADOS)
	g++ -c -g -O0  -std=c++11 test_chain_mod.cpp -o test_chain.o -I$(ACADOS_EXT) -I$(ACADOS) -I$(QPDUNESINCLUDEDIR) -I$(ACADOS_EXT)/eigen
	g++ -g -O0  test_chain.o read_matrix.o read_ocp_qp_in.o model.o jac_chain_nm2.o jac_chain_nm3.o jac_chain_nm4.o  vde_chain_nm2.o vde_chain_nm3.o vde_chain_nm4.o -o test_chain.exe $(LDFLAGS) -lacados -lblasfeo -lqpdunes $(LIBS)
	# valgrind --leak-check=full --show-leak-kinds=all ./test_chain.exe
	./test_chain.exe

read:
	gcc -c -g -O0 ../../test/test_utils/read_ocp_qp_in.c -o read_ocp_qp_in.o -I$(ACADOS) -I$(ACADOS_EXT)
	g++ -c -g -O0 readtst.cpp -o readtst.o -I$(ACADOS)
	g++ -g -O0 read_ocp_qp_in.o readtst.o -o read.exe $(LDFLAGS) -lacados -lblasfeo

simple: simple.c
	gcc -g simple.c -o simple.exe -I$(ACADOS)
	./simple.exe

clean:
	rm *.exe *.o
