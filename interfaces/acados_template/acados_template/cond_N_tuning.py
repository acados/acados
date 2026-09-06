"""Measured choice of `qp_solver_cond_N` for an `AcadosOcpSolver`.

acados condenses an N-stage OCP QP to `qp_solver_cond_N` stages before the
HPIPM solve (PARTIAL_CONDENSING_HPIPM). The default is N (no condensing). The
cost over cond_N is roughly unimodal -- condensing work grows as cond_N falls,
Riccati work grows as it rises -- and the minimum depends on N, nx, nu AND the
machine, so no rule of thumb places it. Measured on acados' chain-of-masses
example (nx 15-63, N 40-100) the default costs 1.7-4.3x in QP time over the
best cond_N; on pure OCP QPs (nx 16-128, N 40-640) 3-13x.

This module chooses cond_N by measuring, through the existing
`AcadosOcpSolver.update_qp_solver_cond_N`, which re-creates the solver at a new
horizon without code generation. Two entry points:

  tune_qp_solver_cond_N(solver, ...)
      Offline: after the first solve, re-solve the SAME problem at candidate
      horizons from a reset iterate (a converged problem re-solved from its own
      iterate does no SQP iteration at all), time the QP per SQP iteration,
      commit the best, restore the iterate. A handful of solves, no closed loop.

  CondNTuner(solver, ...)
      Online: probe on the first solves of a closed loop (`before_solve()` /
      `after_solve()` around each solve), then commit. For loops where an
      offline pass is not possible.

Both use a discrete unimodal search over candidate horizons {1, 2, 4, ..., N}
(about 2*log2(log2 N) probes), take QP time per SQP iteration as the metric
(the only thing cond_N changes; per iteration so that transients do not bias
the probes), and keep the current cond_N unless a candidate is better by more
than `tolerance` (default 10%): below that the difference is timing noise.

No model of the machine, no tuned performance constant.
"""
import math
from typing import Callable, Iterable, List, Optional

__all__ = ['candidates', 'tune_qp_solver_cond_N', 'CondNTuner']


def candidates(N: int) -> List[int]:
    """Powers of two up to N, and N."""
    grid = sorted({2**k for k in range(int(math.log2(N))+1)} | {N})
    return [c for c in grid if 1 <= c <= N]


class _Bracket:
    """Discrete unimodal search: evaluate both ends and the middle of the index
    range, keep the half holding the smaller end, repeat; the middle stays in
    both halves so nothing measured is wasted."""

    def __init__(self, grid: List[int]):
        self.grid = grid
        self.cost = {}
        self.lo, self.hi = 0, len(grid)-1

    def next(self) -> Optional[int]:
        while True:
            lo, hi = self.lo, self.hi
            if hi-lo <= 2:
                for i in range(lo, hi+1):
                    if self.grid[i] not in self.cost:
                        return self.grid[i]
                return None
            mid = (lo+hi)//2
            for i in (lo, hi, mid):
                if self.grid[i] not in self.cost:
                    return self.grid[i]
            if self.cost[self.grid[lo]] <= self.cost[self.grid[hi]]:
                self.hi = mid
            else:
                self.lo = mid

    def report(self, n2: int, seconds: float):
        self.cost[n2] = min(seconds, self.cost.get(n2, float('inf')))

    def best(self) -> int:
        return min(self.cost, key=self.cost.get)


def _require_partial_condensing(solver):
    qp_solver = str(getattr(solver.ocp.solver_options, 'qp_solver', ''))
    if not qp_solver.startswith('PARTIAL_CONDENSING'):
        raise ValueError(f'qp_solver_cond_N can only be changed in place for a partial-condensing '
                         f'QP solver; this solver uses {qp_solver!r}')


def _qp_seconds_per_iteration(solver) -> Optional[float]:
    """QP time per SQP iteration of the last solve; None if no QP was solved
    (a converged problem re-solved from its own iterate does zero iterations)."""
    iters = int(solver.get_stats('sqp_iter'))
    if iters < 1:
        return None
    return float(solver.get_stats('time_qp'))/iters


def _snapshot(solver):
    """The solver's current iterate, to restore after probing (API-version tolerant)."""
    if hasattr(solver, 'store_iterate_to_flat_obj'):
        return ('flat', solver.store_iterate_to_flat_obj())
    if hasattr(solver, 'store_iterate_to_obj'):
        return ('obj', solver.store_iterate_to_obj())
    N = solver.N
    return ('manual', [(solver.get(i, 'x'), solver.get(i, 'u') if i < N else None) for i in range(N+1)])


def _restore(solver, snap):
    kind, it = snap
    if kind == 'flat':
        solver.load_iterate_from_flat_obj(it)
    elif kind == 'obj':
        solver.load_iterate_from_obj(it)
    else:
        for i, (x, u) in enumerate(it):
            solver.set(i, 'x', x)
            if u is not None:
                solver.set(i, 'u', u)


def _cold_start(solver, x0):
    """Reset the iterate, then put every stage at the current initial state.
    A plain reset leaves the trajectory at zero, which for many models is a
    singular configuration (the chain of masses: all masses on one point,
    spring force NaN) and makes the SQP abort before its first QP."""
    solver.reset()
    for i in range(solver.N+1):
        solver.set(i, 'x', x0)


def _probe(solver, n2: int, x0) -> float:
    """One measured solve at horizon n2 from the SAME cold start for every
    candidate, so the SQP has to iterate (a converged problem re-solved from its
    own iterate does no QP at all) and the probes compare."""
    solver.update_qp_solver_cond_N(int(n2))
    _cold_start(solver, x0)
    status = solver.solve()
    value = _qp_seconds_per_iteration(solver)
    if value is None:
        raise RuntimeError(f'tune_qp_solver_cond_N: the solve did no SQP iteration (status {status}); '
                           'is the initial-state constraint set (solve_for_x0 / lbx_0)?')
    return value


def tune_qp_solver_cond_N(solver, repeats: int = 1, tolerance: float = 0.1,
                          grid: Optional[Iterable[int]] = None, verbose: bool = False) -> int:
    """Choose and set `qp_solver_cond_N` by re-solving the current problem.

    Call after at least one `solve()`, with the initial state and guesses in
    place: every probe re-solves the same OCP from the same (converged)
    iterate, so the probes are cheap and comparable. The solver is left at the
    chosen horizon with the same solution as before.

    :param repeats: solves per candidate; the fastest counts (1 is usually enough)
    :param tolerance: keep the current cond_N unless a candidate is faster by
        more than this fraction -- differences below it are timing noise
    :param grid: candidate horizons; default {1, 2, 4, ..., N} plus N
    :returns: the cond_N now set on the solver
    """
    _require_partial_condensing(solver)
    N = solver.N
    current = int(solver.ocp.solver_options.qp_solver_cond_N or N)
    search = _Bracket(list(grid) if grid is not None else candidates(N))
    if verbose:
        print(f'tune_qp_solver_cond_N: N={N}, current cond_N={current}, candidates={search.grid}')
    snap = _snapshot(solver)
    x0 = solver.get(0, 'x')       # stage 0 is pinned to the initial state
    while (n2 := search.next()) is not None:
        for _ in range(repeats):
            search.report(n2, _probe(solver, n2, x0))
        if verbose:
            print(f'  cond_N={n2:4d}: {search.cost[n2]*1e3:8.3f} ms per QP')
    if current not in search.cost:
        search.report(current, _probe(solver, current, x0))
    best = search.best()
    if search.cost[best] < (1.-tolerance)*search.cost[current]:
        choice = best
    else:
        choice = current          # nothing clearly better: do not churn
    solver.update_qp_solver_cond_N(int(choice))
    _restore(solver, snap)        # leave the solver where the caller left it
    if verbose:
        print(f'  -> cond_N={choice} ({search.cost[current]/search.cost[choice]:.2f}x faster QP than cond_N={current})')
    return int(choice)


class CondNTuner:
    """Online tuning inside a closed loop.

        tuner = CondNTuner(solver)
        for each sample:
            tuner.before_solve()
            u = solver.solve_for_x0(x)
            tuner.after_solve()
        tuner.choice   # the committed cond_N once tuner.done

    The probes are real solves whose results the loop uses like any other; the
    metric is QP time per SQP iteration, so the transient does not bias them.
    """

    def __init__(self, solver, tolerance: float = 0.1, grid: Optional[Iterable[int]] = None):
        _require_partial_condensing(solver)
        self.solver = solver
        self.N = solver.N
        self.current = int(solver.ocp.solver_options.qp_solver_cond_N or self.N)
        self.search = _Bracket(list(grid) if grid is not None else candidates(self.N))
        self.tolerance = tolerance
        self.pending: Optional[int] = None
        self.choice: Optional[int] = None
        self.probes = 0

    @property
    def done(self) -> bool:
        return self.choice is not None

    def before_solve(self):
        if self.done:
            return
        n2 = self.search.next()
        if n2 is None:
            self._commit()
            return
        self.solver.update_qp_solver_cond_N(int(n2))
        self.pending = n2

    def after_solve(self):
        if self.pending is None:
            return
        value = _qp_seconds_per_iteration(self.solver)
        if value is None:
            return                # no QP this step (already converged): probe again next step
        self.search.report(self.pending, value)
        self.pending = None
        self.probes += 1
        if self.search.next() is None:
            self._commit()

    def _commit(self):
        best = self.search.best()
        cur = self.search.cost.get(self.current)
        if cur is not None and not self.search.cost[best] < (1.-self.tolerance)*cur:
            best = self.current
        self.choice = int(best)
        self.solver.update_qp_solver_cond_N(self.choice)
