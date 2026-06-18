from anderson_convergence_experiment import create_solver

from acados_template import ACADOS_INFTY
import numpy as np

def timings_experiment(variant, with_abs_cost):
    n_rep = 100
    solver = create_solver(variant, tol=1e-8, with_abs_cost=with_abs_cost)
    initial_guess = solver.get_flat_iterate()

    for anderson_activation_threshold in [0.0, ACADOS_INFTY]:
        time_tot = np.zeros((n_rep,))
        time_glob = np.zeros((n_rep,))
        for i in range(n_rep):
            solver.options_set("anderson_activation_threshold", anderson_activation_threshold)
            solver.set_iterate(initial_guess)
            status = solver.solve()
            time_tot[i] = solver.get_stats("time_tot")
            time_glob[i] = solver.get_stats("time_glob")
        n_iter = solver.get_stats("nlp_iter")
        avg_timing_per_iter = np.min(time_tot) / n_iter * 1e6
        avg_time_glob_per_iter = np.min(time_glob) / n_iter * 1e6
        if anderson_activation_threshold == 0.0:
            time_glob_base = avg_time_glob_per_iter
        else:
            print(f"addtional time for AA per iteration: {avg_time_glob_per_iter - time_glob_base:.4f} microseconds")
        label = "without Anderson" if anderson_activation_threshold == 0 else "with Anderson"
        print(f"{label} - average time per iteration: {avg_timing_per_iter:.4f}, glob: {avg_time_glob_per_iter:.4f} in microseconds over {n_rep} repetitions")

if __name__ == "__main__":
    timings_experiment("EXACT", with_abs_cost=True)
