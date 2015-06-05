CC = gcc
CFLAGS = -O2 -march=native

OBJS = test_hpmpc_blas.o

test_problem: $(OBJS)
	$(CC) -o test.out $(OBJS) -lhpmpc -lm
	@echo
	@echo " Test problem build complete"
	@echo

run: test_problem
	./test.out

clean:
	rm -f *.o
	rm -f test.out

