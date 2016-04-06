SRC = tools.c ocp_qp_solver.c
OBJ = tools.o ocp_qp_solver.o

CFLAGS = -I. -I../hpmpc/include
CFLAGS += -std=c99 -O2 -Wall -Wextra -Werror

all: static_library

static_library:
	$(CC) -c $(SRC) $(CFLAGS)
	ar rcs libacados.a $(OBJ) ../hpmpc/libhpmpc.a

clean:
	rm -f *.o
	rm -f libacados.a
