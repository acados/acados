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

import numpy as np
import casadi
import scipy.interpolate
import matplotlib.pyplot as plt

class TrackSpline:
    def __init__(self, control_points:np.ndarray, spline_seg_length:np.ndarray) -> None:
        if control_points.shape[0] != spline_seg_length.size:
            raise ValueError('control_points and spline_seg_length must have the same length')
        self._num_control_points = control_points.shape[0]
        self._control_points = control_points.copy()
        self._spline_seg_length = spline_seg_length.copy()
        self._break_points = np.cumsum(self._spline_seg_length)
        self._spline_length = self._break_points[-1]
        self._break_points /= self._spline_length
        self._track_spline = scipy.interpolate.CubicSpline(self._break_points, self._control_points)
        self._setup_casadi_functions()

    def _setup_casadi_functions(self):
        spline_coeffs = self._track_spline.c
        self._ca_spline_coeff_x = casadi.MX(spline_coeffs[:,:,0])
        self._ca_spline_coeff_y = casadi.MX(spline_coeffs[:,:,1])
        self._ca_control_points = casadi.MX(self._track_spline.x)
        s = casadi.MX.sym('s')
        seg_idx = casadi.low(self._ca_control_points, s)
        delta_s = s - self._ca_control_points[seg_idx, 0]
        cubic_vector = casadi.vertcat(delta_s**3, delta_s**2, delta_s, 1)
        coeff_x = self._ca_spline_coeff_x[:, seg_idx]
        coeff_y = self._ca_spline_coeff_y[:, seg_idx]
        pos_x = coeff_x.T @ cubic_vector
        pos_y = coeff_y.T @ cubic_vector
        self._ca_fun_pos_x = casadi.Function('pos_x', [s], [pos_x])
        self._ca_fun_pos_y = casadi.Function('pos_y', [s], [pos_y])
        pos_x_1st_deriv = casadi.gradient(pos_x, s)
        pos_y_1st_deriv = casadi.gradient(pos_y, s)
        self._ca_fun_pos_x_1st_deriv = casadi.Function('pos_x_1st_deriv', [s], [pos_x_1st_deriv])
        self._ca_fun_pos_y_1st_deriv = casadi.Function('pos_y_1st_deriv', [s], [pos_y_1st_deriv])
        pos_y_2nd_deriv = casadi.gradient(pos_y_1st_deriv, s)
        pos_x_2nd_deriv = casadi.gradient(pos_x_1st_deriv, s)
        self._ca_fun_pos_x_2nd_deriv = casadi.Function('pos_x_2nd_deriv', [s], [pos_x_2nd_deriv])
        self._ca_fun_pos_y_2nd_deriv = casadi.Function('pos_y_2nd_deriv', [s], [pos_y_2nd_deriv])

        v_s = casadi.MX.sym('v_s')
        a_s = casadi.MX.sym('a_s')
        v_t = v_s * self._spline_length
        a_t = a_s * self._spline_length
        heading = casadi.arctan2(pos_y_1st_deriv, pos_x_1st_deriv)
        omega_t = v_s * casadi.gradient(heading, s)
        alpha_t = a_s * casadi.gradient(heading, s) + v_s * casadi.gradient(omega_t, s)
        self._ca_fun_heading = casadi.Function('heading', [s], [heading])
        self._ca_fun_v_t = casadi.Function('v_t', [s, v_s], [v_t])
        self._ca_fun_a_t = casadi.Function('a_t', [s, v_s, a_s], [a_t])
        self._ca_fun_omega_t = casadi.Function('omega_t', [s, v_s], [omega_t])
        self._ca_fun_alpha_t = casadi.Function('alpha_t', [s, v_s, a_s], [alpha_t])

    @property
    def ca_fun_pos_x(self) -> casadi.Function:
        return self._ca_fun_pos_x

    @property
    def ca_fun_pos_y(self) -> casadi.Function:
        return self._ca_fun_pos_y

    @property
    def ca_fun_pos_x_1st_deriv(self) -> casadi.Function:
        return self._ca_fun_pos_x_1st_deriv

    @property
    def ca_fun_pos_y_1st_deriv(self) -> casadi.Function:
        return self._ca_fun_pos_y_1st_deriv

    @property
    def ca_fun_pos_x_2nd_deriv(self) -> casadi.Function:
        return self._ca_fun_pos_x_2nd_deriv

    @property
    def ca_fun_pos_y_2nd_deriv(self) -> casadi.Function:
        return self._ca_fun_pos_y_2nd_deriv

    @property
    def ca_fun_heading(self) -> casadi.Function:
        return self._ca_fun_heading

    @property
    def ca_fun_v_t(self) -> casadi.Function:
        return self._ca_fun_v_t

    @property
    def ca_fun_a_t(self) -> casadi.Function:
        return self._ca_fun_a_t

    @property
    def ca_fun_omega_t(self) -> casadi.Function:
        return self._ca_fun_omega_t

    @property
    def ca_fun_alpha_t(self) -> casadi.Function:
        return self._ca_fun_alpha_t

    @property
    def spline_length(self) -> float:
        return self._spline_length

    def plot_spline(self):
        s = np.linspace(0, self._break_points[-1], 200)
        scipy_spline = self._track_spline(s, 0)
        casadi_spline = np.zeros((s.size, 2))
        for i in range(s.size):
            casadi_spline[i, 0] = self.ca_fun_pos_x(s[i]).full()
            casadi_spline[i, 1] = self.ca_fun_pos_y(s[i]).full()
        fig = plt.figure(0)
        ax = fig.subplots()
        ax.scatter(self._control_points[:,0], self._control_points[:,1], label='points')
        ax.plot(casadi_spline[:,0], casadi_spline[:,1], label='casadi')
        ax.plot(scipy_spline[:,0], scipy_spline[:,1], label='scipy', linestyle='--')
        for i in range(0, s.size, 5):
            tangnent_vec = np.array([self.ca_fun_pos_x_1st_deriv(s[i]).full()[0, 0],
                                     self.ca_fun_pos_y_1st_deriv(s[i]).full()[0, 0]])
            if np.linalg.norm(tangnent_vec) > 0:
                tangnent_vec /= np.linalg.norm(tangnent_vec)
            ax.arrow(casadi_spline[i,0], casadi_spline[i,1],
                     tangnent_vec[0], tangnent_vec[1])
        ax.legend(loc='upper left')
        ax.set_aspect('equal')

    def plot_trajectories_over_time(self, s:np.ndarray, v_s:np.ndarray, a_s:np.ndarray, timestamps:np.ndarray):
        n_points = timestamps.size
        heading = np.full((n_points,), np.nan)
        v = np.full((n_points,), np.nan)
        a = np.full((n_points,), np.nan)
        omega = np.full((n_points,), np.nan)
        alpha = np.full((n_points,), np.nan)
        for i in range(n_points):
            heading[i] = self.ca_fun_heading(s[i]).full()[0, 0]
            v[i] = self.ca_fun_v_t(s[i], v_s[i]).full()[0, 0]
            omega[i] = self.ca_fun_omega_t(s[i], v_s[i]).full()[0, 0]
            if i < n_points-1:
                a[i] = self.ca_fun_a_t(s[i], v_s[i], a_s[i]).full()[0, 0]
                alpha[i] = self.ca_fun_alpha_t(s[i], v_s[i], a_s[i]).full()[0, 0]
        fig = plt.figure(1)
        ax = fig.subplots(6, 1, sharex=True)
        ax[0].plot(timestamps, s)
        ax[0].set_ylabel('s')
        ax[1].plot(timestamps, heading)
        ax[1].set_ylabel('heading')
        ax[2].plot(timestamps, v)
        ax[2].set_ylabel('v')
        ax[3].plot(timestamps, a)
        ax[3].set_ylabel('a')
        ax[4].plot(timestamps, omega)
        ax[4].set_ylabel('omega')
        ax[5].plot(timestamps, alpha)
        ax[5].set_xlabel('timestamps')
        ax[5].set_ylabel('alpha')

if __name__ == '__main__':
    control_points = np.array([[-2.0, 2.0],
                               [0.0, 0.0],
                               [2.0, 0.0],
                               [4.0, 0.0],
                               [6.0, 0.0],
                               [8.0, 2.0],
                               [8.0, 4.0],
                               [8.0, 6.0],
                               [8.0, 8.0]])
    spline_seg_length = np.array([0.0, np.pi, 2.0, 2.0, 2.0, np.pi, 2.0, 2.0, 2.0])
    track_spline = TrackSpline(control_points, spline_seg_length)
    track_spline.plot_spline()

    v_s_0 = 0.5
    a_s_constant = -0.1
    if not np.isclose(0., a_s_constant):
        tf = (np.sqrt(v_s_0 * v_s_0 + 2 * a_s_constant) - v_s_0) / a_s_constant
    else:
        tf = 1.0 / v_s_0
    n_points = 1000
    timestamps = np.linspace(0, tf, n_points+1, endpoint=True)
    a_s = np.ones(n_points) * a_s_constant
    v_s = v_s_0 + a_s_constant * timestamps
    s = v_s_0*timestamps + 0.5*a_s_constant*timestamps**2
    track_spline.plot_trajectories_over_time(s=s, v_s=v_s, a_s=a_s, timestamps=timestamps)
    plt.show()