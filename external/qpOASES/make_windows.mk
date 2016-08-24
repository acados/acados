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
##	Filename:  make_windows.mk
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.1embedded
##	Date:      2007-2015
##


##
##	definitions for compiling with Visual Studio under Windows
##


################################################################################
# user configuration

# include directories, relative
IDIR =   ${TOP}/include
SRCDIR = ${TOP}/src
BINDIR = ${TOP}/bin

# Matlab include directory (ADAPT TO YOUR LOCAL SETTINGS!)
#MATLAB_IDIR   = ${HOME}/Programs/matlab/extern/include/
MATLAB_IDIR = /usr/local/matlab/extern/include/
MATLAB_LIBDIR = /usr/local/matlab/bin/glnxa64/

# choose between static (1) and dynamic (0) library
MAKE_STATIC_LIB = 1


################################################################################
# do not touch this

CPP = cl
CC  = cl
AR  = ar
RM  = rm
ECHO = echo
CD = cd
CP = copy

# file extensions
OBJEXT = obj
LIBEXT = lib
DLLEXT = so
EXE = .exe
MEXOCTEXT = mex
DEF_TARGET = 
SHARED = /LD

# 32 or 64 depending on target platform
BITS = $(shell getconf LONG_BIT)

# decide on MEX interface extension
ifeq ($(BITS), 32)
	MEXEXT = mexglx
else
	MEXEXT = mexa64
endif

CFLAGS = -nologo -EHsc -DWIN32 -Dsnprintf=_snprintf
#-g -D__DEBUG__ -D__NO_COPYRIGHT__ -D__SUPPRESSANYOUTPUT__

# libraries to link against when building qpOASES_e .so files
LINK_LIBRARIES = 

ifeq ($(MAKE_STATIC_LIB), 1)
	# how to link against the qpOASES_e static library
	QPOASES_LINK = /I${BINDIR} -lqpOASES_e
	# link dependencies when creating executables
	LINK_DEPENDS = ${BINDIR}/libqpOASES_e.${LIBEXT}
else
	# how to link against the qpOASES_e shared library
	QPOASES_LINK = /I${BINDIR} /WL /link ${BINDIR}/libqpOASES_e.lib
	# link dependencies when creating executables
	LINK_DEPENDS = ${BINDIR}/libqpOASES_e.${DLLEXT}
endif



##
##	end of file
##
