include ../../Makefile.rule

OBJS =
OBJS += external_function_interface.o
OBJS += dense_qp_interface.o
OBJS += ocp_nlp_interface.o
OBJS += ocp_qp_interface.o
OBJS += condensing_interface.o
OBJS += options_interface.o
OBJS += sim_interface.o

static_library: $(OBJS)
	ar rcs libacados_c.a $(OBJS)
	@echo
	@echo " libacados_c.a static library build complete."
	@echo

shared_library: $(OBJS)
	gcc -L../../lib -shared -o libacados_c.so $(OBJS) -lacore -lhpipm -lblasfeo -lm
	@echo
	@echo " libacados_c.so shared library build complete."
	@echo

clean:
	rm -f *.o
	rm -f *.s
	rm -f libacados_c.a
	rm -f libacados_c.so
