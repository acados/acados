OBJS = ./src/tools.o ./src/ocp_qp_solver.o

CFLAGS = -O2 -I../../hpmpc/include -I../include -DHPMPC
CFLAGS += -std=c99 -Wall -Wextra -Werror

all: static_library

static_library:
	( cd src ; $(MAKE) obj )
	ar rcs libacados.a $(OBJS) ../hpmpc/libhpmpc.a

clean:
	( cd src ; $(MAKE) clean )
	rm libacados.a
