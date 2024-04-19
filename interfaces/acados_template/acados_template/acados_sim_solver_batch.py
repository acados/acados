from .acados_sim_solver import AcadosSimSolver
from typing import List
from ctypes import (POINTER, c_int, c_void_p)

class AcadosSimSolverBatch():

    def __init__(self, sim_solvers: List[AcadosSimSolver]):

        self.sim_solvers = sim_solvers
        self.N_batch = len(sim_solvers)
        self.shared_lib = sim_solvers[0].shared_lib
        self.model_name = sim_solvers[0].model_name
        self.sim_solvers_pointer = (c_void_p * self.N_batch)()

        for i in range(self.N_batch):
            self.sim_solvers_pointer[i] = self.sim_solvers[i].capsule

        getattr(self.shared_lib, f"{self.model_name}_acados_sim_solve").argtypes = [POINTER(c_void_p), c_int]
        getattr(self.shared_lib, f"{self.model_name}_acados_sim_solve").restype = c_void_p


    def solve(self):
        getattr(self.shared_lib, f"{self.model_name}_acados_sim_batch_solve")(self.sim_solvers_pointer, self.N_batch)

