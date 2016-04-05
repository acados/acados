all: hpmpc/libhpmpc.a

# build HPMPC library
hpmpc/libhpmpc.a:
	( cd hpmpc; $(MAKE) static_library)

# build HPMPC test executable
hpmpc/test_d_ip_hard.o:
	( cd hpmpc; $(MAKE) test_problem)

# build acados test executable
test_problems/test.out:
	( cd test_problems; $(MAKE))

# run the tests
.PHONY: test
test: hpmpc/test_d_ip_hard.o test_problems/test.out
	( cd hpmpc; $(MAKE) run)
	./test_problems/test.out

.PHONY: clean
clean:
	( cd hpmpc; $(MAKE) clean)
	( cd test_problems; $(MAKE) clean)
