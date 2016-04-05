all: hpmpc/libhpmpc.a

# build HPMPC library
hpmpc/libhpmpc.a:
	( cd hpmpc; $(MAKE) static_library)

# build HPMPC test executable
hpmpc/test_d_ip_hard.o:
	( cd hpmpc; $(MAKE) test_problem)

# build acados test executable
acados/test/test.out:
	( cd acados/test; $(MAKE))

# run the tests
.PHONY: test
test: hpmpc/test_d_ip_hard.o acados/test/test.out
	( cd hpmpc; $(MAKE) run)
	./acados/test/test.out


.PHONY: clean
clean:
	( cd hpmpc; $(MAKE) clean)
	( cd acados/test; $(MAKE) clean)


# run a linter
.PHONY: lint
lint: lint-acados lint-hpmpc

ACADOS_LINT_SRC = $(shell find acados -type f -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp')

# TODO: remove these and clean them up
ACADOS_STYLE_FILTER = \
	-build/header_guard, \
	-build/include, \
	-legal/copyright, \
	-readability/casting, \
	-readability/todo, \
	-whitespace/blank_line, \
	-whitespace/braces, \
	-whitespace/comma, \
	-whitespace/comments, \
	-whitespace/end_of_line, \
	-whitespace/line_length, \
	-whitespace/newline, \
	-whitespace/operators, \
	-whitespace/parens, \
	-whitespace/semicolon, \
	-whitespace/tab, \

.PHONY: lint-acados
lint-acados: $(ACADOS_LINT_SRC)
	./cpplint.py --filter="$(ACADOS_STYLE_FILTER)" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 $(ACADOS_LINT_SRC)

# TODO: use CPPLINT.cfg files in subdirectories if hpmpc becomes part of acados
HPMPC_LINT_SRC = $(shell find hpmpc -type f -name '*.c' -o -name '*.cpp' -o -name '*.h' -o -name '*.hpp')

# TODO: remove these and clean them up
HPMPC_STYLE_FILTER = \
	-build/header_guard, \
	-build/include, \
	-build/include_order, \
	-legal/copyright, \
	-readability/braces, \
	-readability/casting, \
	-readability/fn_size, \
	-readability/multiline_comment, \
	-readability/todo, \
	-runtime/int, \
	-whitespace/blank_line, \
	-whitespace/braces, \
	-whitespace/comma, \
	-whitespace/comments, \
	-whitespace/end_of_line, \
	-whitespace/indent, \
	-whitespace/line_length, \
	-whitespace/newline, \
	-whitespace/operators, \
	-whitespace/parens, \
	-whitespace/semicolon, \
	-whitespace/tab, \
	-whitespace/todo, \

.PHONY: lint-hpmpc
lint-hpmpc: $(HPMPC_LINT_SRC)
	./cpplint.py --filter="$(HPMPC_STYLE_FILTER)" --counting=detailed --extensions=c,cpp,h,hpp --linelength=100 $(HPMPC_LINT_SRC)
