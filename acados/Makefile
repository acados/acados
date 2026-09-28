#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# Licensed under the 2-Clause BSD License.


include ../Makefile.rule

obj:
	(cd ocp_nlp; $(MAKE) obj )
	(cd dense_qp; $(MAKE) obj )
	(cd ocp_qp; $(MAKE) obj )
	(cd sim; $(MAKE) obj )
	(cd utils; $(MAKE) obj )

clean:
	(cd ocp_nlp; $(MAKE) clean )
	(cd dense_qp; $(MAKE) clean )
	(cd ocp_qp; $(MAKE) clean )
	(cd sim; $(MAKE) clean )
	(cd utils; $(MAKE) clean )
