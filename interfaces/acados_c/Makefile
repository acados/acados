include ../../Makefile.rule

OBJS =

OBJS += dense_qp.o
OBJS += legacy_create.o
OBJS += ocp_lin.o
# OBJS += ocp_nlp.o
OBJS += ocp_qp.o
OBJS += options.o
OBJS += sim.o

static_library: $(OBJS)
	ar rcs libacados_c.a $(OBJS)
	@echo
	@echo " libacados_c.a static library build complete."
	@echo

clean:
	rm -f *.o
	rm -f *.s
	rm -f libacados_c.a