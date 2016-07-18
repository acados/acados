/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES_e.h
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.1embedded
 *	\date 2007-2015
 */


#if defined(__SINGLE_OBJECT__)

#include <MessageHandling.c>
#include <Utils.c>
#include <Indexlist.c>
#include <Bounds.c>
#include <Constraints.c>
#include <Flipper.c>
#include <Matrices.c>
#include <Options.c>

/*#if (QPOASES_NCMAX == 0)*/
  #include <QProblemB.c>
/*#else*/
  #include <QProblem.c>
/*#endif*/

  /*#include <OQPinterface.c>*/

#else

/*#if (QPOASES_NCMAX == 0)*/
  #include <qpOASES_e/QProblemB.h>
/*#else*/
  #include <qpOASES_e/QProblem.h>
/*#endif*/

  #include <qpOASES_e/extras/OQPinterface.h>

#endif
