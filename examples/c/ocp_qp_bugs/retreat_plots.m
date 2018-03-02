


%% qpOASES warmstart

close all

nmasses = 7;
solver = 'qpOASES_e_N2';
warmstart = 1;

T = plot_comparison(nmasses, solver, warmstart);
T.Interpreter = 'latex';
T.String = 'qpOASES (e) $N^2$ - warmstart';
% exportfig('~/Desktop/ski_retreat/figures/qpOASES_N2.pdf')

% #1: TODO COMPARE WITH COLD START QPOASES TOO?

%% qpOASES coldstart

close all

nmasses = 7;
solver = 'qpOASES_e_N2';
warmstart = 0;

T = plot_comparison(nmasses, solver, warmstart);
T.Interpreter = 'latex';
T.String = 'qpOASES (e) $N^2$ - coldstart';

% exportfig('~/Desktop/ski_retreat/figures/qpOASES_N2_cold.pdf')

%% qpOASES vs dense HPIPM

close all

nmasses = 7;
solver = 'qpOASES_e_N2';
warmstart = 1;

T = plot_comparison(nmasses, solver, warmstart, 'dense_hpipm');
T.Interpreter = 'latex';
T.String = 'qpOASES (e) $N^2$ - warmstart vs hpipm';

% exportfig('~/Desktop/ski_retreat/figures/qpOASES_N2_and_HPIPM.pdf')

%% qpDUNES B0

close all

nmasses = 7;
solver = 'qpDUNES_B0';
warmstart = 1;

T = plot_comparison(nmasses, solver, warmstart);
T.Interpreter = 'latex';
T.String = 'qpDUNES (clipping) - warmstart';

% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B0.pdf')

%% qpDUNES B0 - coldstart

close all

nmasses = 7;
solver = 'qpDUNES_B0';
warmstart = 1;

T = plot_comparison(nmasses, solver, warmstart, 'qpdunes_coldstart');
T.Interpreter = 'latex';
T.String = 'qpDUNES (clipping) - warmstart vs coldstart';

% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B0_cold.pdf')

%% qpDUNES B0 - shifting

close all

nmasses = 7;
solver = 'qpDUNES_B0';
warmstart = 1;

T = plot_comparison(nmasses, solver, warmstart, 'qpdunes_shift');
T.Interpreter = 'latex';
T.String = 'qpDUNES (clipping) - warmstart vs shifting';

% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B0_shift.pdf')

%% qpDUNES B10

close all

nmasses = 7;
solver = 'qpDUNES_B10';
warmstart = 1;
T = plot_comparison(nmasses, solver, warmstart)
T.Interpreter = 'latex';
T.String = 'qpDUNES (B$10$) - warmstart';

% exportfig('~/Desktop/ski_retreat/figures/qpDUNES_B10.pdf')

%% HPMPC B0

close all

nmasses = 7;
solver = 'HPMPC_B0';
warmstart = 0;
T = plot_comparison(nmasses, solver, warmstart);
T.Interpreter = 'latex';
T.String = 'HPMPC - coldstart';

% exportfig('~/Desktop/ski_retreat/figures/HPMPC_B0.pdf')

%% HPMPC B10 (INTERNAL CONDENSING) -- MUST BE WRONG

nmasses = 7;
solver = 'HPMPC_B10';
warmstart = 0;
plot_comparison(nmasses, solver, warmstart)
% exportfig('~/Desktop/ski_retreat/figures/HPMPC_B10_internal_cond.pdf')


%% HPMPC_B10 with acados condensing

close all

nmasses = 7;
solver = 'HPMPC_B10';
warmstart = 0;

T = plot_comparison(nmasses, solver, warmstart, 'hpmpc_ext_cond');
T.Interpreter = 'latex';
T.String = 'qpDUNES (clipping) - warmstart vs coldstart';

