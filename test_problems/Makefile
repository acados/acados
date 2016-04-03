###################################################################################################
#                                                                                                 #
# This file is part of HPMPC.                                                                     #
#                                                                                                 #
# HPMPC -- Library for High-Performance implementation of solvers for MPC.                        #
# Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.                #
#                                                                                                 #
# HPMPC is free software; you can redistribute it and/or                                          #
# modify it under the terms of the GNU Lesser General Public                                      #
# License as published by the Free Software Foundation; either                                    #
# version 2.1 of the License, or (at your option) any later version.                              #
#                                                                                                 #
# HPMPC is distributed in the hope that it will be useful,                                        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                  #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            #
# See the GNU Lesser General Public License for more details.                                     #
#                                                                                                 #
# You should have received a copy of the GNU Lesser General Public                                #
# License along with HPMPC; if not, write to the Free Software                                    #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  #
#                                                                                                 #
# Author: Gianluca Frison, giaf (at) dtu.dk                                                       #
#                                                                                                 #
###################################################################################################

include ./Makefile.rule

ifeq ($(QP_SOLVER), HPMPC)
OBJS = tools.o test_acados_ocp_qp_interface.o
else
OBJS = dummy_main.o
endif

CC = gcc
ifeq ($(QP_SOLVER), HPMPC)
CFLAGS = -O2 -I/opt/hpmpc/include -DHPMPC
else
CFLAGS = -O2 
endif

obj: $(OBJS)
ifeq ($(QP_SOLVER), HPMPC)
	$(CC) -o test.out $(OBJS) -L. /opt/hpmpc/lib/libhpmpc.a -lm
else
	$(CC) -o test.out $(OBJS) -lm
endif

run:
	./test.out

clean:
	rm -f *.o
	rm -f test.out


