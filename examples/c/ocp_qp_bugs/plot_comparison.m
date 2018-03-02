
function [THANDLE, FHANDLE] = plot_comparison(nmasses, solver, warmstart, second_solver, FHANDLE)

if nargin < 4
    second_solver = [];
end

if nargin < 5
    FHANDLE = figure;
end

XLIM = [10 80];

% XXX: horizon length
% YYY: ACADO or acados and cpu_times or iters (or sol_error)
filename_gen = sprintf('ocp_qp_data_nmasses_%d_nsteps_XXX_solver_%s_warmstart_%d_YYY.txt', nmasses, solver, warmstart);

N_values = 10:10:100;

acados_cpu_times = [];
ACADO_cpu_times = [];
second_acados_cpu_times = [];
acados_iters = [];
second_acados_iters = [];
ACADO_iters = [];
errors = [];
second_errors = [];

for ii = 1:length(N_values)
    
    N = N_values(ii);
    
    % load acados timings
    filename_with_N = strrep(filename_gen, 'XXX', num2str(N));

    fhandle_acados_cpu_times = fopen(strrep(filename_with_N, 'YYY', 'acados_cpu_times'), 'r');
    
    acados_cpu_times = [acados_cpu_times; fscanf(fhandle_acados_cpu_times, '%f')'];
 
    fclose(fhandle_acados_cpu_times);

    
    if ~isempty(second_solver)
        % load second acados_timings
        fhandle_second_acados_cpu_times = fopen(strrep(filename_with_N, 'YYY', sprintf('acados_%s_cpu_times', second_solver)), 'r');

        second_acados_cpu_times = [second_acados_cpu_times; fscanf(fhandle_second_acados_cpu_times, '%f')'];

        fclose(fhandle_second_acados_cpu_times);
    end
    
    % load ACADO timings
    fhandle_ACADO_cpu_times = fopen(strrep(filename_with_N, 'YYY', 'ACADO_cpu_times'), 'r');
    
    ACADO_cpu_times = [ACADO_cpu_times; fscanf(fhandle_ACADO_cpu_times, '%f')'];
    
    fclose(fhandle_ACADO_cpu_times);

    % load acados iters
   
    fhandle_acados_iters = fopen(strrep(filename_with_N, 'YYY', 'acados_iters'), 'r');
    
    acados_iters = [acados_iters; fscanf(fhandle_acados_iters, '%f')'];
    
    fclose(fhandle_acados_iters);
    
    % load ACADO iters
    fhandle_ACADO_iters = fopen(strrep(filename_with_N, 'YYY', 'ACADO_iters'), 'r');
    
    ACADO_iters = [ACADO_iters; fscanf(fhandle_ACADO_iters, '%f')'];
    
    fclose(fhandle_ACADO_iters);
    
    if ~isempty(second_solver)
        fhandle_second_acados_iters = fopen(strrep(filename_with_N, 'YYY', sprintf('acados_%s_iters', second_solver)), 'r');

        second_acados_iters = [second_acados_iters; fscanf(fhandle_second_acados_iters, '%f')'];

        fclose(fhandle_second_acados_iters);
    end

    % load errors
    fhandle_errors = fopen(strrep(filename_with_N, 'YYY', 'sol_error'), 'r');
    
    errors = [errors; fscanf(fhandle_errors, '%f')'];
    
    fclose(fhandle_errors);
    
    if ~isempty(second_solver)
        fhandle_second_errors = fopen(strrep(filename_with_N, 'YYY', sprintf('%s_sol_error', second_solver)), 'r');

        second_errors = [second_errors; fscanf(fhandle_second_errors, '%f')'];

        fclose(fhandle_second_errors);
    end
end

acados_cpu_times(isnan(ACADO_cpu_times)) = NaN;

max_ACADO_cpu_times  = nan(length(N_values),1);
max_acados_cpu_times = nan(length(N_values),1);
max_ACADO_iters      = nan(length(N_values),1);
max_acados_iters     = nan(length(N_values),1);
max_error            = nan(length(N_values),1);

if ~isempty(second_solver)
    max_second_acados_cpu_times = nan(length(N_values),1);
    max_second_acados_iters     = nan(length(N_values),1);
    max_second_error            = nan(length(N_values),1);
end

for ii = 1:length(N_values)
    % find worst-time execution time with ACADO
    [max_ACADO_cpu_times(ii), worst_time_pos] = max(ACADO_cpu_times(ii,:));
    % find corresponding time with acados
    max_acados_cpu_times(ii) = acados_cpu_times(ii, worst_time_pos);
    if ~isempty(second_solver)
        max_second_acados_cpu_times(ii) = second_acados_cpu_times(ii, worst_time_pos);
    end
    % find number of iterations
    max_ACADO_iters(ii) = ACADO_iters(ii, worst_time_pos);
    max_acados_iters(ii) = acados_iters(ii, worst_time_pos);
    if ~isempty(second_solver)
        max_second_acados_iters(ii) = second_acados_iters(ii, worst_time_pos);   
    end    
    % find error of same problem instance
    max_error(ii) = errors(ii, worst_time_pos);
    if ~isempty(second_solver)
        max_second_error(ii) = second_errors(ii, worst_time_pos);   
    end
end

yyaxis left

h(1) = subplot(3,1,1)
plot(N_values, 1000*max_ACADO_cpu_times, '-ob', 'linewidth', 1.5)
hold on
plot(N_values, 1000*max_acados_cpu_times, '-or', 'linewidth', 1.5)
if ~isempty(second_solver)
hold on
plot(N_values, 1000*max_second_acados_cpu_times, '-om', 'linewidth', 1.5)
end
ylabel('cpu time [ms]')

set(gca, 'fontsize',18);

yyaxis right

plot(N_values, max_ACADO_cpu_times./max_ACADO_cpu_times, '--b');
hold on
plot(N_values, max_ACADO_cpu_times./max_acados_cpu_times, '--r');
if ~isempty(second_solver)
    hold on
    plot(N_values, max_ACADO_cpu_times./max_second_acados_cpu_times, '--m');
end
xlim(XLIM);

ylabel('speedup')

if ~isempty(second_solver)
    solver_name = second_solver;
    solver_name(solver_name == '_') = ' ';
    legend('ACADO', 'acados', solver_name)
else
    legend('ACADO', 'acados')
end

solver_name = solver;
solver_name(solver_name == '_') = ' ';
THANDLE = title(solver_name);

h(2) = subplot(3,1,2);

% stairs(N_values, max_ACADO_iters, '-ob', 'linewidth', 1.5)
% hold on
% stairs(N_values, max_acados_iters, '-or', 'linewidth', 1.5)
% legend('ACADO', 'acados')
% ylabel('iterations')

% plot(N_values, max_acados_iters - max_ACADO_iters, 'ob', 'linewidth', 1.5)
% hold on
% plot(N_values, zeros(size(max_ACADO_iters)), '--r', 'linewidth', 1.5)

if isempty(second_solver)
    hb = bar(N_values, [max_ACADO_iters max_acados_iters],'grouped');
else
    hb = bar(N_values, [max_ACADO_iters max_acados_iters max_second_acados_iters],'grouped');
end

set(hb(1), 'FaceColor','b')
set(hb(2), 'FaceColor','r')
if ~isempty(second_solver)
    set(hb(3), 'FaceColor','m')
end

xlim([XLIM(1)-5 XLIM(2)+5]);

% xlabel('N');
ylabel('iters')
set(gca, 'fontsize',18);

h(3) = subplot(3,1,3);
semilogy(N_values, max_error, '-ob', 'linewidth', 1.5)
xlim(XLIM);

set(gca, 'fontsize',18);

if ~isempty(second_solver)
    hold on
    semilogy(N_values, max_second_error, '-om', 'linewidth', 1.5)
end
ylabel('error')
xlabel('N');


% [x, y, w, h]
% keyboard

hpos = cell(3);
hpos{1} = get(h(1), 'position');
hpos{2} = get(h(2), 'position');
hpos{3} = get(h(3), 'position');

set(h(1), 'position', hpos{1}.*[1 0.7 1 2]);
set(h(2), 'position', hpos{2}.*[1 0.7 1 0.5]);
set(h(3), 'position', hpos{3}.*[1 1 1 0.5]);

FHANDLE.Position = [100 300 600 500];

end