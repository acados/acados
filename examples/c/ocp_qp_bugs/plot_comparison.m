
function plot_comparison(nmasses, solver, warmstart)

% XXX: horizon length
% YYY: ACADO or acados and cpu_times or iters (or sol_error)
filename_gen = sprintf('ocp_qp_data_nmasses_%d_nsteps_XXX_solver_%s_warmstart_%d_YYY.txt', nmasses, solver, warmstart);

N_values = 10:10:100;

acados_cpu_times = [];
ACADO_cpu_times = [];
acados_iters = [];
ACADO_iters = [];
errors = [];

for ii = 1:length(N_values)
    
    N = N_values(ii);
    
    % load acados timings
    filename_with_N = strrep(filename_gen, 'XXX', num2str(N));

    fhandle_acados_cpu_times = fopen(strrep(filename_with_N, 'YYY', 'acados_cpu_times'), 'r');
    
    acados_cpu_times = [acados_cpu_times; fscanf(fhandle_acados_cpu_times, '%f')'];
 
    fclose(fhandle_acados_cpu_times);

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
    
    % load errors
    fhandle_errors = fopen(strrep(filename_with_N, 'YYY', 'sol_error'), 'r');
    
    errors = [errors; fscanf(fhandle_errors, '%f')'];
    
    fclose(fhandle_errors);
       
end

acados_cpu_times(isnan(ACADO_cpu_times)) = NaN;

max_ACADO_cpu_times  = nan(length(N_values),1);
max_acados_cpu_times = nan(length(N_values),1);
max_ACADO_iters      = nan(length(N_values),1);
max_acados_iters     = nan(length(N_values),1);
max_error            = nan(length(N_values),1);

for ii = 1:length(N_values)
    % find worst-time execution time with ACADO
    [max_ACADO_cpu_times(ii), worst_time_pos] = max(ACADO_cpu_times(ii,:));
    % find corresponding time with acados
    max_acados_cpu_times(ii) = acados_cpu_times(ii, worst_time_pos);
    % find number of iterations
    max_ACADO_iters(ii) = ACADO_iters(ii, worst_time_pos);
    max_acados_iters(ii) = acados_iters(ii, worst_time_pos);
    % find error of same problem instance
    max_error(ii) = errors(ii, worst_time_pos);
end

close all

subplot(3,1,1)
plot(N_values, 1000*max_ACADO_cpu_times, '-ob', 'linewidth', 1.5)
hold on
plot(N_values, 1000*max_acados_cpu_times, '-or', 'linewidth', 1.5)
hold on
plotyy(N_values, nan(size(max_acados_cpu_times)), N_values, [max_ACADO_cpu_times./max_acados_cpu_times ones(size(max_acados_cpu_times))])
xlabel('N');
ylabel('cpu time [ms]')
legend('ACADO', 'acados')
title(solver)

subplot(3,1,2)

% stairs(N_values, max_ACADO_iters, '-ob', 'linewidth', 1.5)
% hold on
% stairs(N_values, max_acados_iters, '-or', 'linewidth', 1.5)
% legend('ACADO', 'acados')
% ylabel('iterations')
stairs(N_values, max_acados_iters - max_ACADO_iters, '-ob', 'linewidth', 1.5)
xlabel('N');
ylabel('iteration error (acados-acado)')

subplot(3,1,3)
semilogy(N_values, max_error, '-ob', 'linewidth', 1.5)
ylabel('error')

end