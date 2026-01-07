#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#


from typing import List, Union
import numpy as np
import casadi as ca
from .utils import cast_to_2d_nparray
from copy import deepcopy
from collections import OrderedDict
from dataclasses import dataclass

@dataclass
class AcadosParam:
    name: str
    value: np.ndarray

class AcadosParamManager:
    """
    Class to manage acados parameters.
    """

    def __init__(self, params: List[AcadosParam], N_horizon: int = 0, use_SX: bool = True):
        """
        Initialize the AcadosParamManager.
        :param params: List of parameters as (name, value) tuples.
        :param N_horizon: Number of stages in the horizon.
        :param use_SX: Whether to use SX or MX for CasADi expressions.
        """

        self._param_values = [OrderedDict()]  # list of dicts for each stage
        self._param_expressions = OrderedDict()

        self._params = params

        if use_SX:
            symbolics = ca.SX.sym
        else:
            symbolics = ca.MX.sym

        for p in params:
            p.value = cast_to_2d_nparray(p.value, p.name)
            self._param_values[0][p.name] = p.value
            self._param_expressions[p.name] = symbolics(p.name, self._param_values[0][p.name].shape)

        # copy default values to all stages
        for _ in range(N_horizon):
            self._param_values.append(deepcopy(self._param_values[0]))

        self._N_horizon = N_horizon


    @property
    def N_horizon(self) -> int:
        return self._N_horizon


    @N_horizon.setter
    def N_horizon(self, N_horizon: int):
        """
        Set the horizon length and adjust parameter values accordingly.
        """
        if not isinstance(N_horizon, int) or N_horizon < 0:
            raise ValueError("N_horizon must be a non-negative integer")

        # If increasing, add new stages with copies stage 0
        old_N = len(self._param_values) - 1
        if N_horizon > old_N:
            for _ in range(N_horizon - old_N):
                self._param_values.append(deepcopy(self._param_values[0]))

        # Initialize with default values
        for p in self._params:
            for n in range(self._N_horizon, N_horizon):
                self._param_values[n][p.name] = p.value

        self._N_horizon = N_horizon


    def get_value(self, name: str, stage: int) -> np.ndarray:
        """
        Get the value of a parameter by name.

        :param name: Name of the parameter.
        :param stage: Stage index.
        :return: Value of the parameter as a numpy array.
        """
        if stage > self._N_horizon or stage < 0:
            raise IndexError(f"Stage index {stage} out of bounds for horizon length {self._N_horizon}.")
        return self._param_values[stage][name]


    def get_expression(self, name: str) -> str:
        """
        Get the expression of a parameter by name.

        :param name: Name of the parameter.
        :return: CasADi expression of the parameter.
        """
        return self._param_expressions[name]


    def set_value(self, name: str, value: np.ndarray, stage: int):
        """
        Set the value of a parameter by name.

        :param name: Name of the parameter.
        :param value: New value of the parameter as a numpy array.
        """
        if stage > self._N_horizon or stage < 0:
            raise IndexError(f"Stage index {stage} out of bounds for horizon length {self._N_horizon}.")
        self._param_values[stage][name] = value


    def get_p_stagewise_values(self, stage: int) -> np.ndarray:
        """
        Get the values of the stage-wise parameters for a given stage.

        :param stage: stage index.
        :return: numpy array of the parameter vector for a given stage.
        """
        if stage > self._N_horizon or stage < 0:
            raise IndexError(f"Stage index {stage} out of bounds for horizon length {self._N_horizon}.")
        return ca.vertcat(*[ca.vec(v) for v in self._param_values[stage].values()]).full()


    def get_p_stagewise_values_flat(self,) -> np.ndarray:
        """
        Get all values of the stage-wise parameter as a flat vector.

        :return: numpy array of the parameter vector for a given stage.
        """
        return ca.vertcat(*[self.get_p_stagewise_values(stage) for stage in range(self._N_horizon + 1)]).full()


    def get_p_stagewise_expression(self) -> Union[ca.SX, ca.MX]:
        """
        Get the stage-wise expression of a parameter by name.

        :return: CasADi expression of the parameter vector.
        """
        return ca.vertcat(*[ca.vec(v) for v in self._param_expressions.values()])
