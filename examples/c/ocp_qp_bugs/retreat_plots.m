


%% qpOASES 

nmasses = 7;
solver = 'qpOASES_e_N2';
warmstart = 1;

plot_comparison(nmasses, solver, warmstart)
% exportfig('~/Desktop/ski_retreat/figures/qpoases_N2.pdf')

% TODO COMPARE WITH COLD START QPOASES TOO?

%% qpOASES vs dense HPIPM

nmasses = 7;
solver = 'qpOASES_e_N2';
warmstart = 1;

plot_comparison(nmasses, solver, warmstart, 'dense_hpipm')
% exportfig('~/Desktop/ski_retreat/figures/qpoases_N2_and_hpipm.pdf')

%% qpDUNES B0

nmasses = 7;
solver = 'qpDUNES_B0';
warmstart = 1;

plot_comparison(nmasses, solver, warmstart)
% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B0.pdf')

%% qpDUNES B10

nmasses = 7;
solver = 'qpDUNES_B10';
warmstart = 1;
plot_comparison(nmasses, solver, warmstart)
% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B10.pdf')

%% HPMPC B0

nmasses = 7;
solver = 'HPMPC_B0';
warmstart = 0;
plot_comparison(nmasses, solver, warmstart)
% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B10.pdf')