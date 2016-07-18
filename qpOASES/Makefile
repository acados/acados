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
##	Filename:  Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.1embedded
##	Date:      2007-2015
##

include make.mk

##
##	targets
##


all: src examples

src:
	@cd $@; ${MAKE} -s 

examples: src
	@cd $@; ${MAKE} -s

doc:
	@cd $@; ${MAKE} -s 

testing: src
	@cd testing/c; ${MAKE} -s

clean:
	@cd src               && ${MAKE} -s clean
	@cd examples          && ${MAKE} -s clean
	@cd bin               && ${RM} -f *.* *{EXE}
	@cd testing/c         && ${MAKE} -s clean

clobber: clean

.PHONY : all src examples doc testing clean clobber


##
##   end of file
##
