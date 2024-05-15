from .acados_sim_solver import AcadosSimSolver
from .acados_sim import AcadosSim
from typing import List
from ctypes import (POINTER, c_int, c_void_p)


class AcadosSimBatchSolver():
    """
    Batch Integrator for parallel integration.

        :param sim: type :py:class:`~acados_template.acados_sim.AcadosSim`
        :param N_batch: batch size, positive integer
        :param json_file: Default: 'acados_sim.json'
        :verbose: bool, default: True
    """

    __sim_solvers : List[AcadosSimSolver]

    def __init__(self, sim: AcadosSim, N_batch: int, json_file: str = 'acados_sim.json', verbose: bool=True):

        if not isinstance(N_batch, int) or N_batch <= 0:
            raise Exception("AcadosSimBatchSolver: argument N_batch should be a positive integer.")

        self.__N_batch = N_batch
        self.__sim_solvers = [AcadosSimSolver(sim, json_file=json_file, build=n==0, generate=n==0, verbose=verbose) for n in range(self.N_batch)]

        self.__shared_lib = self.sim_solvers[0].shared_lib
        self.__model_name = self.sim_solvers[0].model_name
        self.__sim_solvers_pointer = (c_void_p * self.N_batch)()

        for i in range(self.N_batch):
            self.__sim_solvers_pointer[i] = self.sim_solvers[i].capsule

        getattr(self.__shared_lib, f"{self.__model_name}_acados_sim_batch_solve").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.__shared_lib, f"{self.__model_name}_acados_sim_batch_solve").restype = c_void_p

        if not self.sim_solvers[0].acados_lib_uses_omp:
            print("Warning: Please compile the acados shared library with openmp and the number of threads set to 1, i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1.")


    def solve(self):
        """
        Solve the simulation problem with current input for all `N_batch` integrators.
        """
        getattr(self.__shared_lib, f"{self.__model_name}_acados_sim_batch_solve")(self.__sim_solvers_pointer, self.__N_batch)


    @property
    def sim_solvers(self):
        """List of AcadosSimSolvers."""
        return self.__sim_solvers

    @property
    def N_batch(self):
        """Batch size."""
        return self.__N_batch

