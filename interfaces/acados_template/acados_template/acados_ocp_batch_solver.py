from .acados_ocp_solver import AcadosOcpSolver
from .acados_ocp import AcadosOcp
from typing import List
from ctypes import (POINTER, c_int, c_void_p)


class AcadosOcpBatchSolver():
    """
    Batch Integrator for parallel integration.

        :param sim: type :py:class:`~acados_template.acados_sim.AcadosOcp`
        :param N_batch: batch size, positive integer
        :param json_file: Default: 'acados_sim.json'
        :verbose: bool, default: True
    """

    __ocp_solvers : List[AcadosOcpSolver]

    def __init__(self, ocp: AcadosOcp, N_batch: int, json_file: str = 'acados_ocp.json', verbose: bool=True):

        if not isinstance(N_batch, int) or N_batch <= 0:
            raise Exception("AcadosOcpBatchSolver: argument N_batch should be a positive integer.")

        self.__N_batch = N_batch
        self.__ocp_solvers = [AcadosOcpSolver(ocp, json_file=json_file, build=n==0, generate=n==0, verbose=verbose) for n in range(self.N_batch)]

        self.__shared_lib = self.ocp_solvers[0].shared_lib
        self.__name = self.ocp_solvers[0].name
        self.__ocp_solvers_pointer = (c_void_p * self.N_batch)()

        for i in range(self.N_batch):
            self.__ocp_solvers_pointer[i] = self.ocp_solvers[i].capsule

        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve").restype = c_void_p

        if self.ocp_solvers[0].acados_lib_uses_omp:
            msg = "Note: Please make sure that the acados shared library is compiled with the number of threads set to 1,\n"
        else:
            msg = "Warning: Please compile the acados shared library with openmp and the number of threads set to 1,\n"

        msg += "i.e. with the flags -DACADOS_WITH_OPENMP=ON -DACADOS_NUM_THREADS=1.\n" + \
                   "See https://github.com/acados/acados/pull/1089 for more details."
        print(msg)


    def solve(self):
        """
        Call solve for all `N_batch` solvers.
        """
        getattr(self.__shared_lib, f"{self.__name}_acados_batch_solve")(self.__ocp_solvers_pointer, self.__N_batch)


    @property
    def ocp_solvers(self):
        """List of AcadosOcpSolvers."""
        return self.__ocp_solvers

    @property
    def N_batch(self):
        """Batch size."""
        return self.__N_batch

