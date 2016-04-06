SRC = tools.c ocp_qp_solver.c

CFLAGS = -O2 -I../../hpmpc/include -I../include -DHPMPC
CFLAGS += -std=c99 -Wall -Wextra -Werror

all: obj

obj: $(SRC)
	($(CC) -c $(SRC) $(CFLAGS) )

clean:
	rm -f *.o
