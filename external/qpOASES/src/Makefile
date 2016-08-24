##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  src/Makefile
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.1embedded
##	Date:      2007-2015
##

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR}

QPOASES_OBJECTS = \
	QProblem.${OBJEXT} \
	QProblemB.${OBJEXT} \
	Bounds.${OBJEXT} \
	Constraints.${OBJEXT} \
	Indexlist.${OBJEXT} \
	Flipper.${OBJEXT} \
	Utils.${OBJEXT} \
	Options.${OBJEXT} \
	Matrices.${OBJEXT} \
	MessageHandling.${OBJEXT}

QPOASES_EXTRAS_OBJECTS = \
	OQPinterface.${OBJEXT}

QPOASES_DEPENDS = \
	${IDIR}/qpOASES_e.h \
	${IDIR}/qpOASES_e/QProblem.h \
	${IDIR}/qpOASES_e/QProblemB.h \
	${IDIR}/qpOASES_e/Flipper.h \
	${IDIR}/qpOASES_e/Bounds.h \
	${IDIR}/qpOASES_e/Constraints.h \
	${IDIR}/qpOASES_e/Indexlist.h \
	${IDIR}/qpOASES_e/Utils.h \
	${IDIR}/qpOASES_e/Constants.h \
	${IDIR}/qpOASES_e/Types.h \
	${IDIR}/qpOASES_e/Options.h \
	${IDIR}/qpOASES_e/Matrices.h \
	${IDIR}/qpOASES_e/MessageHandling.h \
	${IDIR}/qpOASES_e/UnitTesting.h


##
##	targets
##

all: ${LINK_DEPENDS}


${BINDIR}/libqpOASES_e.${LIBEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating static lib" $@
	@${AR} r $@ $^

${BINDIR}/libqpOASES_e.${DLLEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating shared lib" $@
	@${CC} ${DEF_TARGET} ${SHARED} $^ ${LINK_LIBRARIES}

clean:
	@${ECHO} "Cleaning up (src)"
	@${RM} -f *.${OBJEXT} *.${LIBEXT} *.${DLLEXT}

clobber: clean


%.${OBJEXT}: %.c ${QPOASES_DEPENDS}
	@${ECHO} "Creating" $@
	${CC} ${DEF_TARGET} -c ${IFLAGS} ${CFLAGS} $<


##
##	end of file
##
