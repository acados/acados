function opts = set_integrator_opts(q, stepsize, n_steps, forw_sens_mode, adj_sens_mode,butcher_mode)
    if nargin < 6
        butcher_mode = 'irk';
    end
    if butcher_mode == 'irk'
        disp(['using ' num2str(q) ' Gauss-Legendre nodes per step'])
        [~, ~, A_butcher, b, c ] = collo_poly_butcher( q, stepsize );
    else
        disp('using ERK4')
        A_butcher = zeros(4);
        A_butcher(2,1) = .5; A_butcher(3,2) = .5; A_butcher(4,3) = 1;
        b = [1/6, 1/3, 1/3, 1/6]';
        c = [0; .5; .5; 1];
    end
        
    opts.A_butcher = A_butcher;
    opts.b_butcher = b;
    opts.c_butcher = c;
    opts.n_steps = n_steps;

    opts.dt = stepsize/opts.n_steps;
    opts.A_dt = A_butcher * opts.dt;
    opts.b_dt = opts.b_butcher * opts.dt;
    opts.n_stages = length(b);
    opts.forw_sens_mode = forw_sens_mode;
    opts.adj_sens_mode = adj_sens_mode;
end