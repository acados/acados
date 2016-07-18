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
%Returns a struct containing values for all options to be used within qpOASES_e.
%
%Call
%    options = qpOASES_e_options( 'default' );
%    options = qpOASES_e_options( 'reliable' );
%    options = qpOASES_e_options( 'MPC' );
%to obtain a set of default options or a pre-defined set of options tuned
%for reliable or fast QP solution, respectively.
%
%Call
%    options = qpOASES_e_options( 'option1',value1,'option2',value2,... )
%to obtain a set of default options but with 'option1' set to value1 etc.
%
%Call
%    options = qpOASES_e_options( oldOptions,'option1',value1,... )
%to obtain a copy of the options struct oldOptions but with 'option1' set
%to value1 etc.
%
%Call
%    options = qpOASES_e_options( 'default', 'option1',value1,... )
%    options = qpOASES_e_options( 'reliable','option1',value1,... )
%    options = qpOASES_e_options( 'MPC',     'option1',value1,... )
%to obtain a set of default options or a pre-defined set of options tuned
%for reliable or fast QP solution, respectively, but with 'option1' set to 
%value1 etc.
%
%
%qpOASES_e features the following options:
%  maxIter                    -  Maximum number of iterations (if set
%                                to -1, a value is chosen heuristically)
%  maxCpuTime                 -  Maximum CPU time in seconds (if set
%                                to -1, only iteration limit is used)
%  printLevel                 -  0: no printed output,
%                                1: only error messages are printed,
%                                2: iterations and error messages are printed,
%                                3: all available messages are printed.
%
%  enableRamping              -  Enables (1) or disables (0) ramping.
%  enableFarBounds            -  Enables (1) or disables (0) the use of 
%                                far bounds.
%  enableFlippingBounds       -  Enables (1) or disables (0) the use of 
%                                flipping bounds.
%  enableRegularisation       -  Enables (1) or disables (0) automatic 
%                                Hessian regularisation.
%  enableFullLITests          -  Enables (1) or disables (0) condition-hardened 
%                                (but more expensive) LI test.
%  enableNZCTests             -  Enables (1) or disables (0) nonzero curvature 
%                                tests.
%  enableDriftCorrection      -  Specifies the frequency of drift corrections:
%                                0: turns them off, 
%                                1: uses them at each iteration etc.
%  enableCholeskyRefactorisation - Specifies the frequency of a full re-
%                                factorisation of projected Hessian matrix:
%                                0: turns them off, 
%                                1: uses them at each iteration etc.
%  enableEqualities           -  Specifies whether equalities should be treated 
%                                as always active (1) or not (0)
%
%  terminationTolerance       -  Relative termination tolerance to stop homotopy.
%  boundTolerance             -  If upper and lower bounds differ less than this
%                                tolerance, they are regarded equal, i.e. as 
%                                equality constraint.
%  boundRelaxation            -  Initial relaxation of bounds to start homotopy 
%                                and initial value for far bounds.
%  epsNum                     -  Numerator tolerance for ratio tests.
%  epsDen                     -  Denominator tolerance for ratio tests.
%  maxPrimalJump              -  Maximum allowed jump in primal variables in 
%                                nonzero curvature tests.
%  maxDualJump                -  Maximum allowed jump in dual variables in 
%                                linear independence tests.
%
%  initialRamping             -  Start value for ramping strategy.
%  finalRamping               -  Final value for ramping strategy.
%  initialFarBounds           -  Initial size for far bounds.
%  growFarBounds              -  Factor to grow far bounds.
%  initialStatusBounds        -  Initial status of bounds at first iteration:
%                                 0: all bounds inactive,
%                                -1: all bounds active at their lower bound,
%                                +1: all bounds active at their upper bound.
%  epsFlipping                -  Tolerance of squared Cholesky diagonal factor 
%                                which triggers flipping bound.
%  numRegularisationSteps     -  Maximum number of successive regularisation steps.
%  epsRegularisation          -  Scaling factor of identity matrix used for 
%                                Hessian regularisation.
%  numRefinementSteps         -  Maximum number of iterative refinement steps.
%  epsIterRef                 -  Early termination tolerance for iterative 
%                                refinement.
%  epsLITests                 -  Tolerance for linear independence tests.
%  epsNZCTests                -  Tolerance for nonzero curvature tests.
%
%
%See also QPOASES, QPOASES_SEQUENCE, QPOASES_AUXINPUT
%
%
%For additional information see the qpOASES User's Manual or
%visit http://www.qpOASES.org/.
%
%Please send remarks and questions to support@qpOASES.org!
function [ options ] = qpOASES_e_options( varargin )

	firstIsStructOrScheme = 0;

	if ( nargin == 0 ) 
		options = qpOASES_e_default_options();
	else
		if ( isstruct( varargin{1} ) )
			if ( mod( nargin,2 ) ~= 1 )
				error('ERROR (qpOASES_e_options): Options must be specified in pairs!');
			end
			options = varargin{1};
			firstIsStructOrScheme = 1;
		else
			if ( ischar( varargin{1} ) )
				if ( mod( nargin,2 ) == 0 )
					options = qpOASES_e_default_options();
				else
					if ( ( nargin > 1 ) && ( ischar( varargin{nargin} ) ) )
						error('ERROR (qpOASES_e_options): Options must be specified in pairs!');
					end

					switch ( varargin{1} )
						case 'default'
							options = qpOASES_e_default_options();
						case 'reliable'
							options = qpOASES_e_reliable_options();
						case {'MPC','mpc','fast'}
							options = qpOASES_e_MPC_options();
						otherwise
							error( ['ERROR (qpOASES_e_options): Only the following option schemes are defined: ''default'', ''reliable'', ''MPC''!'] );
							
					end
					firstIsStructOrScheme = 1;
				end
			else
				error('ERROR (qpOASES_e_options): First argument needs to be a string or an options struct!');
			end
		end
	end

	% set options to user-defined values  
	for i=(1+firstIsStructOrScheme):2:nargin

		argName  = varargin{i};
		argValue = varargin{i+1};

		if ( ( isempty( argName ) ) || ( ~ischar( argName ) ) )
			error('ERROR (qpOASES_e_options): Argmument no. %d has to be a non-empty string!',i );
        end
			
		if ( ( ischar(argValue) ) || ( ~isscalar( argValue ) ) )
			error('ERROR (qpOASES_e_options): Argmument no. %d has to be a scalar constant!',i+1 );
        end

		if ( ~isfield( options,argName ) )
			error('ERROR (qpOASES_e_options): Argmument no. %d is an invalid option!',i );
		end

		eval( ['options.',argName,' = ',num2str(argValue),';'] );

	end

end


function [ options ] = qpOASES_e_default_options( )

	% setup options struct with default values
	options = struct(	'maxIter',                       -1, ...
						'maxCpuTime',                    -1, ...
						'printLevel',                     1, ...
						...
						'enableRamping',                  1, ...
						'enableFarBounds',                1, ...
						'enableFlippingBounds',           1, ...
						'enableRegularisation',           0, ...
						'enableFullLITests',              0, ...
						'enableNZCTests',                 1, ...
						'enableDriftCorrection',          1, ...
						'enableCholeskyRefactorisation',  0, ...
						'enableEqualities',               0, ...
						...
						'terminationTolerance',           5.0e6*eps, ...
						'boundTolerance',                 1.0e6*eps, ...
						'boundRelaxation',                1.0e4, ...
						'epsNum',                        -1.0e3*eps, ...
						'epsDen',                         1.0e3*eps, ...
						'maxPrimalJump',                  1.0e8, ...
						'maxDualJump',                    1.0e8, ...
						...
    					'initialRamping',                 0.5, ...
						'finalRamping',                   1.0, ...
						'initialFarBounds',               1.0e6, ...
						'growFarBounds',                  1.0e3, ...
						'initialStatusBounds',            -1, ...
						'epsFlipping',                    1.0e3*eps, ...
						'numRegularisationSteps',         0, ...
						'epsRegularisation',              1.0e3*eps, ...
						'numRefinementSteps',             1, ...
						'epsIterRef',                     1.0e2*eps, ...
						'epsLITests',                     1.0e5*eps, ...
						'epsNZCTests',                    3.1e3*eps );

end



function [ options ] = qpOASES_e_reliable_options( )

	% setup options struct with values for most reliable QP solution
	options = qpOASES_e_default_options( );

	options.enableFullLITests             =  1;
	options.enableCholeskyRefactorisation =  1;

	options.numRefinementSteps            =  2;

end


function [ options ] = qpOASES_e_MPC_options( )

	% setup options struct with values for most reliable QP solution
	options = qpOASES_e_default_options( );

	options.enableRamping                 =  0;
	options.enableFarBounds               =  1;
	options.enableFlippingBounds          =  0;
	options.enableRegularisation          =  1;
	options.enableNZCTests                =  0;
	options.enableDriftCorrection         =  0;
	options.enableEqualities              =  1;

	options.terminationTolerance          =  1.0e9*eps;
	
	options.initialStatusBounds           =  0;
	options.numRegularisationSteps        =  1;
	options.numRefinementSteps            =  0;

end
