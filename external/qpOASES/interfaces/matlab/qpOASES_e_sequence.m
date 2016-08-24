%qpOASES -- An Implementation of the Online Active Set Strategy.
%Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
%Christian Kirches et al. All rights reserved.
%
%qpOASES is distributed under the terms of the
%GNU Lesser General Public License 2.1 in the hope that it will be
%useful, but WITHOUT ANY WARRANTY; without even the implied warranty
%of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%See the GNU Lesser General Public License for more details.
%
%---------------------------------------------------------------------------------
%
%qpOASES_e_sequence is intended to solve a sequence of quadratic
%programming (QP) problems of the following form:
%
%                min   1/2*x'Hx + x'g
%                s.t.  lb  <=  x <= ub
%                      lbA <= Ax <= ubA  {optional}
%
%I) Call
%
%    [QP,x,fval,exitflag,iter,lambda,auxOutput] = ...
%               qpOASES_e_sequence( 'i',H,g,A,lb,ub,lbA,ubA{,options{,auxInput}} )
%or 
%    [QP,x,fval,exitflag,iter,lambda,auxOutput] = ...
%               qpOASES_e_sequence( 'i',H,g,lb,ub{,options{,auxInput}} )
%
%for initialising and solving the first above-mentioned QP of the sequence
%starting from an initial guess x0. H must be a symmetric (possibly indefinite) 
%matrix and all vectors g, lb, ub, lbA, ubA have to be given as column vectors. 
%Options can be generated using the qpOASES_e_options command, otherwise default
%values are used. Optionally, further auxiliary inputs may be generated 
%using qpOASES_e_auxInput command and passed to the solver.
%Both matrices H or A may be passed in sparse matrix format.
%
%II) Call
%
%     [x,fval,exitflag,iter,lambda,auxOutput] = ...
%                      qpOASES_e_sequence( 'h',QP,g,lb,ub,lbA,ubA{,options} )
%or
%     [x,fval,exitflag,iter,lambda,auxOutput] = ...
%                      qpOASES_e_sequence( 'h',QP,g,lb,ub{,options} )
%
%for hotstarting from the previous QP solution to the one of the next QP
%given by the vectors g, lb, ub, lbA, ubA. Options can be generated using the 
%qpOASES_e_options command, otherwise default values are used.
%
%III) Call
%
%     [x,fval,exitflag,iter,lambda,auxOutput] = ...
%                     qpOASES_e_sequence( 'm',QP,H,g,A,lb,ub,lbA,ubA{,options} )
%
%for hotstarting from the previous QP solution to the one of the next QP
%given by the matrices H, A and the vectors g, lb, ub, lbA, ubA. The previous
%active set serves as a starting guess. If the new projected Hessian matrix
%turns out to be not positive definite, qpOASES recedes to a safe initial active
%set guess automatically. This can result in a high number of iterations iter.
%Options can be generated using the qpOASES_e_options command, otherwise default
%values are used.
%
%IV) Call
%
%     [x,lambda,workingSetB,workingSetC] = ...
%                     qpOASES_e_sequence( 'e',QP,g,lb,ub,lbA,ubA{,options} )
%
%for solving the equality constrained QP with constraints determined by the
%current active set. All inequalities and bounds which were not active in the
%previous solution might be violated. This command does not alter the internal
%state of qpOASES_e. Instead of calling this command multiple times, it is
%possible to supply several columns simultaneously in g, lb, ub, lbA, and ubA.
%Options can be generated using the qpOASES_e_options command, otherwise default
%values are used.
%
%V) Having solved the last QP of your sequence, call
%
%     qpOASES_e_sequence( 'c',QP )
%
%in order to cleanup the internal memory.
%
%
%Optional outputs (only x is mandatory):
%    x            -  Optimal primal solution vector (if exitflag==0).
%    fval         -  Optimal objective function value (if exitflag==0).
%    exitflag     -   0: QP solved,
%                     1: QP could not be solved within given number of iterations,
%                    -1: QP could not be solved due to an internal error,
%                    -2: QP is infeasible (and thus could not be solved),
%                    -3: QP is unbounded (and thus could not be solved).
%    iter         -  Number of active set iterations actually performed.
%    lambda       -  Optimal dual solution vector (if exitflag==0).
%    auxOutput    -  Struct containing auxiliary outputs as described below.
%
%The auxOutput struct contains the following entries:
%    workingSetB  -  Working set of bounds at point x.
%    workingSetC  -  Working set of constraints at point x.
%                    The working set is a subset of the active set (indices
%                    of bounds/constraints that hold with equality) yielding
%                    a set linearly independent of bounds/constraints.
%                    The working sets are encoded as follows:
%                     1: bound/constraint at its upper bound
%                     0: bound/constraint not at any bound
%                    -1: bound/constraint at its lower bound
%    cpuTime      -  Internally measured CPU time for solving QP problem.
%
%See also QPOASES_E_OPTIONS, QPOASES_E_AUXINPUT, QPOASES_E
%
%
%For additional information see the qpOASES User's Manual or
%visit http://www.qpOASES.org/.
%
%Please send remarks and questions to support@qpOASES.org!
