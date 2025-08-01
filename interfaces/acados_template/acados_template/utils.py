# -*- coding: future_fstrings -*-
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

from typing import Union
import json
import os
import shutil
import sys
import platform
import urllib.request
from subprocess import DEVNULL, STDOUT, call
if os.name == 'nt':
    from ctypes import wintypes
    from ctypes import WinDLL as DllLoader
else:
    from ctypes import CDLL as DllLoader
import numpy as np
from casadi import DM, MX, SX, CasadiMeta, Function
import casadi as ca
from contextlib import contextmanager


TERA_VERSION = "0.2.0"

PLATFORM2TERA = {
    "linux": "linux",
    "darwin": "osx",
    "win32": "windows"
}

ACADOS_INFTY = 1e10

@contextmanager
def set_directory(path: str):
    """Sets the cwd within the context"""
    origin = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def check_if_square(mat: np.ndarray, name: str):
    if mat.shape[0] != mat.shape[1]:
        raise ValueError(f"Matrix {name} must be square, got shape {mat.shape}.")
    return

def get_acados_path():
    ACADOS_PATH = os.environ.get('ACADOS_SOURCE_DIR')
    if not ACADOS_PATH:
        acados_template_path = os.path.dirname(os.path.abspath(__file__))
        acados_path = os.path.join(acados_template_path, '..','..','..')
        ACADOS_PATH = os.path.realpath(acados_path)
        msg = 'Warning: Did not find environment variable ACADOS_SOURCE_DIR, '
        msg += 'guessed ACADOS_PATH to be {}.\n'.format(ACADOS_PATH)
        msg += 'Please export ACADOS_SOURCE_DIR to avoid this warning.'
        print(msg)
    return ACADOS_PATH


def get_python_interface_path():
    ACADOS_PYTHON_INTERFACE_PATH = os.environ.get('ACADOS_PYTHON_INTERFACE_PATH')
    if not ACADOS_PYTHON_INTERFACE_PATH:
        acados_path = get_acados_path()
        ACADOS_PYTHON_INTERFACE_PATH = os.path.join(acados_path, 'interfaces', 'acados_template', 'acados_template')
    return ACADOS_PYTHON_INTERFACE_PATH


def get_tera_exec_path():
    TERA_PATH = os.environ.get('TERA_PATH')
    if not TERA_PATH:
        TERA_PATH = os.path.join(get_acados_path(), 'bin', 't_renderer') + get_binary_ext()

    # convert to absolute path
    TERA_PATH = os.path.abspath(TERA_PATH)
    return TERA_PATH


def acados_lib_is_compiled_with_openmp(acados_lib: DllLoader, verbose: bool) -> bool:
    # find out if acados was compiled with OpenMP
    try:
        acados_lib_uses_omp = getattr(acados_lib, 'omp_get_thread_num') is not None
    except AttributeError as e:
        acados_lib_uses_omp = False
    if verbose:
        if acados_lib_uses_omp:
            print('acados was compiled with OpenMP.')
        else:
            print('acados was compiled without OpenMP.')
    return acados_lib_uses_omp


def get_shared_lib(shared_lib_name: str, winmode = None) -> DllLoader:
    if winmode is not None:
        shared_lib = DllLoader(shared_lib_name, winmode=winmode)
    else:
        # for compatibility with older python versions
        shared_lib = DllLoader(shared_lib_name)
    return shared_lib


def check_casadi_version():
    casadi_version = CasadiMeta.version()
    major_minor = casadi_version.split('.')
    major = int(major_minor[0])
    minor = int(major_minor[1])
    if major < 3 or (major == 3 and minor < 4): # < 3.4
        raise Exception(f'CasADi version {casadi_version} is not supported. '
                        'Please use a version >= 3.4.0.')

    if major > 3 or (major == 3 and minor > 7): # >= 3.7
        print(f"Warning: CasADi version {casadi_version} is not tested with acados yet.")
    elif major == 3 and minor < 7:
        print(f"Warning: Full featured acados requires CasADi version >= 3.7, got {casadi_version}.")


def check_casadi_version_supports_p_global():
    try:
        from casadi import extract_parametric, cse
    except ImportError:
        raise ImportError("CasADi version does not support extract_parametric or cse functions.\nNeeds nightly-se2 release or later, see: https://github.com/casadi/casadi/releases/tag/nightly-se2")


def get_simulink_default_opts() -> dict:
    python_interface_path = get_python_interface_path()
    abs_path = os.path.join(python_interface_path, 'simulink_default_opts.json')
    with open(abs_path , 'r') as f:
        simulink_default_opts = json.load(f)
    return simulink_default_opts


def is_casadi_SX(x):
    if isinstance(x, ca.SX):
        return True
    return False


def is_column(x):
    if isinstance(x, np.ndarray):
        if x.ndim == 1:
            return True
        elif x.ndim == 2 and x.shape[1] == 1:
            return True
        else:
            return False
    elif isinstance(x, (MX, SX, DM)):
        if x.shape[1] == 1:
            return True
        elif x.shape[0] == 0 and x.shape[1] == 0:
            return True
        else:
            return False
    elif x == None or x == []:
        return False
    else:
        raise TypeError("is_column expects one of the following types: np.ndarray, casadi.MX, casadi.SX."
                        + " Got: " + str(type(x)))


def is_empty(x):
    if isinstance(x, (MX, SX, DM)):
        return x.is_empty()
    elif isinstance(x, np.ndarray):
        return True if np.prod(x.shape) == 0 else False
    elif x is None:
        return True
    elif isinstance(x, (set, list, str)):
        return True if len(x) == 0 else False
    elif isinstance(x, (float, int)):
        return False
    else:
        raise TypeError("is_empty expects one of the following types: casadi.MX, casadi.SX, "
                        + "None, numpy array empty list, set. Got: " + str(type(x)))


def casadi_length(x):
    if isinstance(x, (MX, SX, DM)):
        return int(np.prod(x.shape))
    elif x is None:
        return 0
    elif isinstance(x, list):
        return len(x)
    else:
        raise TypeError("casadi_length expects one of the following types: casadi.MX, casadi.SX."
                        + " Got: " + str(type(x)))


def get_shared_lib_ext():
    if sys.platform == 'darwin':
        return '.dylib'
    elif os.name == 'nt':
        return '.dll'
    else:
        return '.so'

def get_shared_lib_dir():
    if os.name == 'nt':
        return 'bin'
    else:
        return 'lib'

def get_shared_lib_prefix():
    if os.name == 'nt':
        return ''
    else:
        return 'lib'

def get_binary_ext():
    if os.name == 'nt':
        return '.exe'
    else:
        return ''

def get_architecture_amd64_arm64():
    # common uname -m results
    # https://en.wikipedia.org/wiki/Uname
    current_arch = platform.machine()
    amd64_compatible = ["i3", "i6", "amd", "x86"]
    arm64_compatible = ["arm", "aarch"]
    if any([current_arch.lower().startswith(arch) for arch in amd64_compatible]):
        return "amd64"
    elif any([current_arch.lower().startswith(arch) for arch in arm64_compatible]):
        return "arm64"
    else:
        raise RuntimeError(f"Your detected architecture {current_arch} may not be compatible with amd64 or arm64.")

def get_tera() -> str:
    tera_path = get_tera_exec_path()
    acados_path = get_acados_path()

    # check if tera exists and is executable
    if os.path.exists(tera_path) and os.access(tera_path, os.X_OK):
        return tera_path

    try:
        arch = get_architecture_amd64_arm64()
    except RuntimeError as e:
        print(e)
        print("Try building tera_renderer from source at https://github.com/acados/tera_renderer")
        sys.exit(1)

    binary_ext = get_binary_ext()
    repo_url = "https://github.com/acados/tera_renderer/releases"
    url = "{}/download/v{}/t_renderer-v{}-{}-{}{}".format(
        repo_url, TERA_VERSION, TERA_VERSION, PLATFORM2TERA[sys.platform], arch, binary_ext)

    manual_install = 'For manual installation follow these instructions:\n'
    manual_install += '1 Download binaries from {}\n'.format(url)
    manual_install += '2 Copy them in {}/bin\n'.format(acados_path)
    manual_install += '3 Strip the version and platform and architecture from the binaries: '
    manual_install += f'as t_renderer-v{TERA_VERSION}-P-A{binary_ext} -> t_renderer{binary_ext})\n'
    manual_install += '4 Enable execution privilege on the file "t_renderer" with:\n'
    manual_install += '"chmod +x {}"\n\n'.format(tera_path)

    msg = "\n"
    msg += 'Tera template render executable not found, '
    msg += 'while looking in path:\n{}\n'.format(tera_path)
    msg += 'In order to be able to render the templates, '
    msg += 'you need to download the tera renderer binaries from:\n'
    msg += '{}\n\n'.format(repo_url)
    msg += 'Do you wish to set up Tera renderer automatically?\n'
    msg += 'y/N? (press y to download tera or any key for manual installation)\n'

    if input(msg) != 'y':
        msg_cancel = "\nYou cancelled automatic download.\n\n"
        msg_cancel += manual_install
        msg_cancel += "Once installed re-run your script.\n\n"
        print(msg_cancel)

        sys.exit(1)

    # check if parent directory exists otherwise create it
    tera_dir = os.path.split(tera_path)[0]
    if not os.path.exists(tera_dir):
        print(f"Creating directory {tera_dir}")
        os.makedirs(tera_dir)

    # Download tera
    print(f"Downloading {url}")
    with urllib.request.urlopen(url) as response, open(tera_path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    print("Successfully downloaded t_renderer.")
    # make executable
    os.chmod(tera_path, 0o755)
    print("Successfully made t_renderer executable.")
    return tera_path


def render_template(in_file, out_file, output_dir, json_path, template_glob=None):

    acados_path = os.path.dirname(os.path.abspath(__file__))
    if template_glob is None:
        head, in_file = os.path.split(in_file)
        template_glob = os.path.join(acados_path, 'c_templates_tera', head, '**', '*')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with set_directory(output_dir):
        tera_path = get_tera()

        # call tera as system cmd
        os_cmd = f"{tera_path} '{template_glob}' '{in_file}' '{json_path}' '{out_file}'"
        # Windows cmd.exe can not cope with '...', so use "..." instead:
        if os.name == 'nt':
            os_cmd = os_cmd.replace('\'', '\"')

        status = os.system(os_cmd)
        if status != 0:
            raise RuntimeError(f'Rendering of {in_file} failed!\n\nAttempted to execute OS command:\n{os_cmd}\n\n')



def casadi_expr_to_string(expr) -> str:
    string = ''
    for ii in range(casadi_length(expr)):
        string += f"{expr[ii,:]}\n"
    return string

## Conversion functions
def make_object_json_dumpable(input):
    if isinstance(input, (np.ndarray)):
        return input.tolist()
    elif isinstance(input, (SX)):
        try:
            return input.serialize()
            # for more readable json output:
            # return casadi_expr_to_string(input)
        except: # for older CasADi versions
            return ''
    elif isinstance(input, (MX)):
        # NOTE: MX expressions can not be serialized, only Functions.
        return input.__str__()
    elif isinstance(input, (DM)):
        return input.full()
    else:
        raise TypeError(f"Cannot make input of type {type(input)} dumpable.")


def format_class_dict(d):
    """
    removes the __ artifact from class to dict conversion
    """
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = format_class_dict(v)

        out_key = k.split('__', 1)[-1]
        out[k.replace(k, out_key)] = v
    return out

def get_default_simulink_opts() -> dict:
    print("get_default_simulink_opts is deprecated, use get_simulink_default_opts instead."
          + " This function will be removed in a future release.")
    return get_simulink_default_opts()


def J_to_idx(J):
    if not isinstance(J, np.ndarray):
        raise TypeError('J_to_idx: J must be a numpy array.')
    if J.ndim != 2:
        raise ValueError('J_to_idx: J must be a 2D numpy array.')
    nrows = J.shape[0]
    idx = np.zeros((nrows, ))
    for i in range(nrows):
        this_idx = np.nonzero(J[i,:])[0]
        if len(this_idx) != 1:
            raise ValueError('Invalid J matrix structure detected, ' \
                'must contain exactly one nonzero element per row.')
        if this_idx.size > 0 and J[i,this_idx[0]] != 1:
            raise ValueError('J matrices can only contain 1 and 0 entries.')
        idx[i] = this_idx[0]
    return idx


def J_to_idx_slack(J):
    nrows = J.shape[0]
    ncol = J.shape[1]
    idx = np.zeros((ncol, ))
    i_idx = 0
    for i in range(nrows):
        this_idx = np.nonzero(J[i,:])[0]
        if len(this_idx) == 1:
            idx[i_idx] = i
            i_idx = i_idx + 1
        elif len(this_idx) > 1:
            raise ValueError('J_to_idx_slack: Invalid J matrix. ' \
                'Found more than one nonzero in row ' + str(i))
        if this_idx.size > 0 and J[i,this_idx[0]] != 1:
            raise ValueError('J_to_idx_slack: J matrices can only contain 1s, ' \
                 'got J(' + str(i) + ', ' + str(this_idx[0]) + ') = ' + str(J[i,this_idx[0]]) )
    if not i_idx == ncol:
            raise ValueError('J_to_idx_slack: J must contain a 1 in every column!')
    return idx


def check_if_nparray_and_flatten(val, name) -> np.ndarray:
    if not isinstance(val, np.ndarray):
        raise TypeError(f"{name} must be a numpy array, got {type(val)}")
    return val.reshape(-1)

def check_if_nparray_or_casadi_symbolic_and_flatten(val, name) -> np.ndarray:
    if not isinstance(val, (np.ndarray, SX, MX)):
        raise Exception(f"{name} must be array of type np.ndarray, casadi.SX, or casadi.MX, got {type(val)}")

    if isinstance(val, (SX, MX)):
        return ca.reshape(val, val.numel(), 1)
    else:
        return val.reshape(-1)


def check_if_2d_nparray(val, name) -> None:
    if not isinstance(val, np.ndarray):
        raise TypeError(f"{name} must be a numpy array, got {type(val)}")
    if val.ndim != 2:
        raise ValueError(f"{name} must be a 2D numpy array, got shape {val.shape}")
    return


def check_if_2d_nparray_or_casadi_symbolic(val, name) -> None:
    if isinstance(val, (SX, MX, DM)):
        return
    if not isinstance(val, np.ndarray):
        raise Exception(f"{name} must be a array of type np.ndarray, casadi.SX, or casadi.MX, got {type(val)}")
    if val.ndim != 2:
        raise Exception(f"{name} must be a 2D array of type np.ndarray, casadi.SX, or casadi.MX, got shape {val.shape}")


def print_J_to_idx_note():
    print("NOTE: J* matrix is converted to zero based vector idx* vector, which is returned here.")


def acados_dae_model_json_dump(model):

    # load model
    x = model.x
    xdot = model.xdot
    u = model.u
    z = model.z
    p = model.p

    f_impl = model.f_impl_expr
    model_name = model.name

    # create struct with impl_dae_fun, casadi_version
    fun_name = model_name + '_impl_dae_fun'
    impl_dae_fun = Function(fun_name, [x, xdot, u, z, p], [f_impl])

    casadi_version = CasadiMeta.version()
    str_impl_dae_fun = impl_dae_fun.serialize()

    dae_dict = {"str_impl_dae_fun": str_impl_dae_fun, "casadi_version": casadi_version}

    # dump
    json_file = model_name + '_acados_dae.json'
    with open(json_file, 'w') as f:
        json.dump(dae_dict, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
    print("dumped ", model_name, " dae to file:", json_file, "\n")


def set_up_imported_gnsf_model(acados_ocp):

    gnsf = acados_ocp.gnsf_model

    # load model
    phi_fun = Function.deserialize(gnsf['phi_fun'])
    phi_fun_jac_y = Function.deserialize(gnsf['phi_fun_jac_y'])
    phi_jac_y_uhat = Function.deserialize(gnsf['phi_jac_y_uhat'])
    get_matrices_fun = Function.deserialize(gnsf['get_matrices_fun'])

    # obtain gnsf dimensions
    size_gnsf_A = get_matrices_fun.size_out(0)
    acados_ocp.dims.gnsf_nx1 = size_gnsf_A[1]
    acados_ocp.dims.gnsf_nz1 = size_gnsf_A[0] - size_gnsf_A[1]
    acados_ocp.dims.gnsf_nuhat = max(phi_fun.size_in(1))
    acados_ocp.dims.gnsf_ny = max(phi_fun.size_in(0))
    acados_ocp.dims.gnsf_nout = max(phi_fun.size_out(0))

    # save gnsf functions in model
    acados_ocp.model.phi_fun = phi_fun
    acados_ocp.model.phi_fun_jac_y = phi_fun_jac_y
    acados_ocp.model.phi_jac_y_uhat = phi_jac_y_uhat
    acados_ocp.model.get_matrices_fun = get_matrices_fun

    # get_matrices_fun = Function([model_name,'_gnsf_get_matrices_fun'], {dummy},...
    #  {A, B, C, E, L_x, L_xdot, L_z, L_u, A_LO, c, E_LO, B_LO,...
    #   nontrivial_f_LO, purely_linear, ipiv_x, ipiv_z, c_LO});
    get_matrices_out = get_matrices_fun(0)
    acados_ocp.model.gnsf_nontrivial_f_LO = int(get_matrices_out[12])
    acados_ocp.model.gnsf_purely_linear = int(get_matrices_out[13])

    if "f_lo_fun_jac_x1k1uz" in gnsf:
        f_lo_fun_jac_x1k1uz = Function.deserialize(gnsf['f_lo_fun_jac_x1k1uz'])
        acados_ocp.model.f_lo_fun_jac_x1k1uz = f_lo_fun_jac_x1k1uz
    else:
        dummy_var_x1 = SX.sym('dummy_var_x1', acados_ocp.dims.gnsf_nx1)
        dummy_var_x1dot = SX.sym('dummy_var_x1dot', acados_ocp.dims.gnsf_nx1)
        dummy_var_z1 = SX.sym('dummy_var_z1', acados_ocp.dims.gnsf_nz1)
        dummy_var_u = SX.sym('dummy_var_z1', acados_ocp.dims.nu)
        dummy_var_p = SX.sym('dummy_var_z1', acados_ocp.dims.np)
        empty_var = SX.sym('empty_var', 0, 0)

        empty_fun = Function('empty_fun', \
            [dummy_var_x1, dummy_var_x1dot, dummy_var_z1, dummy_var_u, dummy_var_p],
                [empty_var])
        acados_ocp.model.f_lo_fun_jac_x1k1uz = empty_fun

    del acados_ocp.gnsf_model


def idx_perm_to_ipiv(idx_perm):
    n = len(idx_perm)
    vec = list(range(n))
    ipiv = np.zeros(n)

    print(n, idx_perm)
    # import pdb; pdb.set_trace()
    for ii in range(n):
        idx0 = idx_perm[ii]
        for jj in range(ii,n):
            if vec[jj]==idx0:
                idx1 = jj
                break
        tmp = vec[ii]
        vec[ii] = vec[idx1]
        vec[idx1] = tmp
        ipiv[ii] = idx1

    ipiv = ipiv-1 # C 0-based indexing
    return ipiv


def print_casadi_expression(f: Union[MX, SX, DM]):
    for ii in range(casadi_length(f)):
        print(f[ii,:])


def verbose_system_call(cmd, verbose=True, shell=False):
    return call(
        cmd,
        stdout=None if verbose else DEVNULL,
        stderr=None if verbose else STDOUT,
        shell=shell
    )
