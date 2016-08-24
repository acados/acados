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
%qpOASES_e solves (a series) of quadratic programming (QP) problems of the
%following form:
%
%                min   1/2*x'Hx + x'g
%                s.t.  lb  <=  x <= ub
%                      lbA <= Ax <= ubA  {optional}
%
%Call
%
%    [x,fval,exitflag,iter,lambda,auxOutput] = 
%                     qpOASES_e( H,g,A,lb,ub,lbA,ubA{,options{,auxInput}} )
%
%for solving the above-mentioned QP. H must be a symmetric (but possibly 
%indefinite) matrix and all vectors g, lb, ub, lbA, ubA have to be given
%as column vectors. Options can be generated using the qpOASES_e_options command, 
%otherwise default values are used. Optionally, further auxiliary inputs
%may be generated using qpOASES_e_auxInput command and passed to the solver.
%Both matrices H or A may be passed in sparse matrix format.
%
%Call
%
%    [x,fval,exitflag,iter,lambda,auxOutput] =
%                     qpOASES_e( H,g,lb,ub{,options{,auxInput}} )
%
%for solving the above-mentioned QP without general constraints.
%
%
%Optional outputs (only x is mandatory):
%    x            -  Optimal primal solution vector (if exitflag==0).
%    fval         -  Optimal objective function value (if exitflag==0).
%    exitflag     -   0: QP problem solved,
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
%
%If not a single QP but a sequence of QPs with varying vectors is to be solved,
%the i-th QP is given by the i-th columns of the QP vectors g, lb, ub, lbA, ubA
%(i.e. they are matrices in this case). Both matrices H and A remain constant.
%
%See also QPOASES_E_OPTIONS, QPOASES_E_AUXINPUT, QPOASES_E_SEQUENCE
%
%
%For additional information see the qpOASES User's Manual or
%visit http://www.qpOASES.org/.
%
%Please send remarks and questions to support@qpOASES.org!
