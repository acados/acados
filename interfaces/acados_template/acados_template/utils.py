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

from typing import Union, Optional
import json
import os
import shutil
import sys
import hashlib
import platform
import urllib.request
import warnings
import subprocess
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

TERA_DEFAULT_VERSION = "0.2.1"

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
        warnings.warn(f"CasADi version {casadi_version} is not tested with acados yet.")
    elif major == 3 and minor < 7:
        warnings.warn(f"Full featured acados requires CasADi version >= 3.7, got {casadi_version}.")


def check_casadi_version_supports_p_global():
    try:
        from casadi import extract_parametric, cse
    except ImportError:
        raise ImportError("CasADi version does not support extract_parametric or cse functions.\nPlease use CasADi >= 3.7.2")

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

def is_none_or_empty_list(x):
    if x is None:
        return True
    elif isinstance(x, list) and len(x) == 0:
        return True
    else:
        return False


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

def is_scalar_integer(x):
    """Check if x is a scalar integer (int or numpy integer type)"""
    return np.isscalar(x) and isinstance(x, (int, np.integer))


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


def get_os_str():
    if sys.platform == 'darwin':
        return 'mac'
    elif os.name == 'nt':
        return 'pc'
    else:
        return 'unix'

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

def _version_tuple(version_str: str):
    """Convert a version string like '0.2.0' to a tuple of ints for comparison."""
    return tuple(int(x) for x in version_str.strip().split('.'))


def is_tera_version_sufficient(tera_path: str, required_version: str) -> bool:
    """Return True if the t_renderer at tera_path meets the required_version.

    Args:
        tera_path: Absolute path to the t_renderer executable.
        required_version: Minimum acceptable version string, e.g. '0.2.0'.

    Old versions (<v0.2.0) do not support the --version flag, which is treated
    as an outdated installation that must be updated.
    """
    try:
        result = subprocess.run(
            [tera_path, '--version'],
            capture_output=True, text=True, timeout=10
        )
    except Exception:
        return False

    if result.returncode != 0:
        # Versions that do not support --version are too old
        return False

    # The output is expected to be e.g. "t_renderer 0.2.0\n" or just "0.2.0\n"
    output = result.stdout.strip()
    # Extract the last whitespace-separated token as the version number
    version_str = output.split()[-1] if output else ""
    try:
        return _version_tuple(version_str) >= _version_tuple(required_version)
    except (ValueError, IndexError):
        warnings.warn(
            f"Could not parse t_renderer version string: '{version_str}'. "
            "Treating the installed version as insufficient and re-downloading.",
            RuntimeWarning
        )
        return False


def get_tera(tera_version: Optional[str] = None, force_download = None) -> str:
    if force_download is not None:
        warnings.warn("get_tera: force_download is deprecated in 0.5.6. Do not pass the option. Now, this function checks if it is available, otherwise attempts downloading or raises an error.", DeprecationWarning, stacklevel=2)

    if tera_version is None:
        tera_version = TERA_DEFAULT_VERSION
    tera_path = get_tera_exec_path()
    acados_path = get_acados_path()

    # check if tera exists, is executable, and has a sufficient version
    if os.path.exists(tera_path) and os.access(tera_path, os.X_OK):
        if is_tera_version_sufficient(tera_path, tera_version):
            return tera_path
        else:
            print(f"\nt_renderer found at {tera_path} but its version does not meet the "
                    f"minimum requirement (>= v{tera_version}). Updating automatically...")

    try:
        arch = get_architecture_amd64_arm64()
    except RuntimeError as e:
        raise RuntimeError(
            f"{e}\nTry building tera_renderer from source at https://github.com/acados/tera_renderer"
        ) from e

    binary_ext = get_binary_ext()
    repo_url = "https://github.com/acados/tera_renderer/releases"
    url = "{}/download/v{}/t_renderer-v{}-{}-{}{}".format(
        repo_url, tera_version, tera_version, PLATFORM2TERA[sys.platform], arch, binary_ext)

    manual_install = 'For manual installation follow these instructions:\n'
    manual_install += '1 Download binaries from {}\n'.format(url)
    manual_install += '2 Copy them in {}/bin\n'.format(acados_path)
    manual_install += '3 Strip the version and platform and architecture from the binaries: '
    manual_install += f'as t_renderer-v{tera_version}-P-A{binary_ext} -> t_renderer{binary_ext})\n'
    manual_install += '4 Enable execution privilege on the file "t_renderer" with:\n'
    manual_install += '"chmod +x {}"\n\n'.format(tera_path)

    try:
        print("Setting up tera renderer")
        # check if parent directory exists otherwise create it
        tera_dir = os.path.split(tera_path)[0]
        if not os.path.exists(tera_dir):
            print(f"Creating directory {tera_dir}")
            os.makedirs(tera_dir)

        # Download tera
        print(f"Downloading {url}")
        with urllib.request.urlopen(url, timeout=60) as response, open(tera_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print("Successfully downloaded t_renderer.")
        # make executable
        os.chmod(tera_path, 0o755)
        print("Successfully made t_renderer executable.")
    except Exception as e:
        msg = "\n"
        msg += 'Tera renderer executable not found, '
        msg += 'while looking in path:\n{}\n'.format(tera_path)
        msg += 'In order to be able to render the templates, '
        msg += 'you need to download the tera renderer binaries from:\n'
        msg += '{}\n\n'.format(repo_url)
        msg += '\nAutomatic installation failed with: {}\n\n'.format(e)
        msg += manual_install
        raise RuntimeError(msg) from e

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
            print(f"\nRendering file {in_file} failed.\n\n",
                  "Known issues:\n",
                  "1) older Linux versions with default tera binaries have issues where a compatible libc.so is not found.\n",
                  "To fix this. Run the following in a Python file:\n" \
                  "from acados_template import get_tera\n",
                  "get_tera()\n\n",
                  "2) ROS templates are not compatibile with old tera version. Only relevant if you try generating a ROS node.")
            raise RuntimeError(f'Rendering of {in_file} failed!\n\nAttempted to execute OS command:\n{os_cmd}\n\n')



def casadi_expr_to_string(expr) -> str:
    string = ''
    for ii in range(casadi_length(expr)):
        string += f"{expr[ii,:]}\n"
    return string

def make_object_json_dumpable(input):
    '''
    Convert numpy arrays and CasADi DM objects to lists before JSON dump.
    NOTE: Serialization of CasADi MX and SX objects requires a StringSerializer and is now implemented in AcadosModel.
    '''
    if isinstance(input, (np.ndarray)):
        return input.tolist()
    elif isinstance(input, (DM)):
        return input.full().tolist()
    elif is_scalar_integer(input):
        return int(input)
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


def J_to_idx(J):
    J = cast_to_2d_nparray(J, 'J')
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
    J = cast_to_2d_nparray(J, 'J')
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


def ns_from_idxs_rev(idxs_rev) -> int:
    if is_empty(idxs_rev):
        return 0
    else:
        ns = int(np.max(idxs_rev) + 1)
        for i in range(ns):
            if i not in idxs_rev:
                raise ValueError(f"Detected ns = {ns}, but i is not in idxs_rev = {idxs_rev}, the slack with index {i} is thus not contained in the problem.")
        return ns

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


def cast_to_1d_nparray(val, name) -> np.ndarray:
    try:
        val = np.asarray(val)
    except:
        raise TypeError(f"Failed to cast {name} to np.array, expected array-like type got {type(val)}.")

    val = np.atleast_1d(np.squeeze(val))

    if val.ndim > 1:
        raise ValueError(f"Expected vector-like array for {name}, got {val.shape}.")

    return val

def use_int_or_cast_to_1d_nparray(val, name) -> Union[int, np.ndarray]:
    if isinstance(val, int):
        return val
    else:
        try:
            return cast_to_1d_nparray(val, name)
        except:
            raise TypeError(f"{name} must be an integer or array-like type, got {type(val)}.")


def cast_to_1d_nparray_or_casadi_symbolic(val, name) -> np.ndarray:
    if isinstance(val, (SX, MX, DM)):
        if val.shape[0] == 1 or val.shape[1] == 1:
            return val
        else:
            raise ValueError(f"Expected vector-like array for {name}, got {val.shape}.")
    else:
        return cast_to_1d_nparray(val, name)


def cast_to_2d_nparray(val, name) -> np.ndarray:
    try:
        val = np.asarray(val)
    except:
        raise TypeError(f"Failed to cast {name} to np.array, expected array-like type got {type(val)}.")

    if val.ndim != 2:
        raise ValueError(f"Expected two dimensional array for {name}, got {val.shape}.")

    return val


def cast_to_2d_nparray_or_casadi_symbolic(val, name) -> np.ndarray:
    if isinstance(val, (SX, MX, DM)):
        return val
    else:
        return cast_to_2d_nparray(val, name)


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

    # create struct with impl_dae_fun, casadi_version
    fun_name = model.name + '_impl_dae_fun'
    impl_dae_fun = Function(fun_name, [x, xdot, u, z, p], [f_impl])

    casadi_version = CasadiMeta.version()
    str_impl_dae_fun = impl_dae_fun.serialize()

    dae_dict = {"str_impl_dae_fun": str_impl_dae_fun, "casadi_version": casadi_version}

    # dump
    json_file = model.name + '_acados_dae.json'
    with open(json_file, 'w') as f:
        json.dump(dae_dict, f, default=make_object_json_dumpable, indent=4, sort_keys=True)
    print("dumped ", model.name, " dae to file:", json_file, "\n")


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

def status_to_str(status):
    status_dict = {
        -1: "ACADOS_UNKNOWN",
        0: "ACADOS_SUCCESS",
        1: "ACADOS_NAN_DETECTED",
        2: "ACADOS_MAXITER",
        3: "ACADOS_MINSTEP",
        4: "ACADOS_QP_FAILURE",
        5: "ACADOS_READY",
        6: "ACADOS_UNBOUNDED",
        7: "ACADOS_TIMEOUT",
        8: "ACADOS_QPSCALING_BOUNDS_NOT_SATISFIED",
        9: "ACADOS_INFEASIBLE"
    }
    return status_dict.get(status, "UNKNOWN_STATUS")

def str_to_status_ipopt(status_str):
    str_dict = {
        "Solve_Succeeded": 0,
        "Solved_To_Acceptable_Level": 0,
        "Maximum_Iterations_Exceeded": 2,
        "Search_Direction_Becomes_Too_Small": 3,
        "Diverging_Iterates": 6,
        "Infeasible_Problem_Detected": 9,
    }
    return str_dict.get(status_str, -1)

OCP_COMPARE_IGNORED_FIELD_PATHS = [
    ('external_function_files_model',),
    ('external_function_files_ocp',),
    ('json_loaded',),
    ('dims', 'n_global_data'),
]

def hash_class_instance(obj) -> str:
    """Create a hash of a class instance based on its attributes."""

    if obj is None:
        class_dict = {}
    elif isinstance(obj, dict):
        class_dict = obj
    elif hasattr(obj, 'to_dict'):
        class_dict = obj.to_dict()
    else:
        class_dict = obj.__dict__

    global OCP_COMPARE_IGNORED_FIELD_PATHS
    for field_path in OCP_COMPARE_IGNORED_FIELD_PATHS:
        child = class_dict
        *path, field_to_remove = field_path
        for p in path:
            child = child.get(p)
            if child is None:
                break
        else:
            child.pop(field_to_remove, None)

    json_str = json.dumps(class_dict, default=make_object_json_dumpable, sort_keys=True)
    hash_md5 = hashlib.md5(json_str.encode('utf-8')).hexdigest()
    # print(f"MD5 hash of the object: {hash_md5}")

    return hash_md5



def compare_ocp_formulations(ocp_1, ocp_2, tol_code_reuse):
    """
    Compare every entry of two OCP objects, ignoring certain fields.

    Args:
        ocp_1: first OCP object with a to_dict() method
        ocp_2: second OCP object with a to_dict() method
        tol_code_reuse: absolute tolerance used when comparing numpy arrays

    Returns:
        List of field paths that do not match
    """
    ocp_1.make_consistent()
    ocp_2.make_consistent()
    dict_1 = ocp_1.to_dict()
    dict_2 = ocp_2.to_dict()

    global OCP_COMPARE_IGNORED_FIELD_PATHS
    for field_path in OCP_COMPARE_IGNORED_FIELD_PATHS:
        _remove_field_path(dict_1, field_path)
        _remove_field_path(dict_2, field_path)

    mismatched_fields = []
    _compare_recursive(dict_1, dict_2, tol_code_reuse, mismatched_fields)

    return mismatched_fields


def _remove_field_path(data, field_path):
    """Remove a nested field given by field_path (a sequence of keys) from data, in place."""
    child = data
    *path, field_to_remove = field_path
    for p in path:
        child = child.get(p)
        if child is None:
            return
    child.pop(field_to_remove, None)


def _compare_recursive(data_1, data_2, tol_code_reuse, mismatched_fields, path=""):
    """
    Recursively compare data_1 and data_2.
    Appends mismatched field paths to mismatched_fields.
    """
    if isinstance(data_1, dict) and isinstance(data_2, dict):
        all_keys = set(data_1) | set(data_2)
        for key in all_keys:
            current_path = f"{path}.{key}" if path else key
            if key not in data_1 or key not in data_2:
                mismatched_fields.append(current_path)
            else:
                _compare_recursive(data_1[key], data_2[key], tol_code_reuse, mismatched_fields, current_path)

    elif isinstance(data_1, (list, tuple)) and isinstance(data_2, (list, tuple)):
        if len(data_1) != len(data_2):
            mismatched_fields.append(path)
        else:
            for i, (item_1, item_2) in enumerate(zip(data_1, data_2)):
                current_path = f"{path}[{i}]"
                _compare_recursive(item_1, item_2, tol_code_reuse, mismatched_fields, current_path)

    else:
        # numpy arrays and CasADi DM objects for comparison
        try:
            # Compare numeric arrays with tolerance
            if isinstance(data_1, (np.ndarray, DM)) or isinstance(data_2, (np.ndarray, DM)):
                # cast DM
                arr_1 = data_1.full() if isinstance(data_1, DM) else data_1
                arr_2 = data_2.full() if isinstance(data_2, DM) else data_2
                if arr_1.shape != arr_2.shape or not np.allclose(arr_1, arr_2, atol=tol_code_reuse, rtol=0.0):
                    mismatched_fields.append(path)
                return
            elif data_1 != data_2:
                mismatched_fields.append(path)

        except (TypeError, ValueError):
            if data_1 != data_2:
                mismatched_fields.append(path)

def verify_weighting_matrix(A, name, tol=1e-10):
    """
    Check if A is square, symmetric, and (positive semidefinite and diagonal) or positive definite matrix.
    Raises an exception otherwise.
    Args:
    A: matrix to check
    name: name of the matrix (for error messages)
    tol: numerical tolerance for symmetry and positive definiteness checks
    """
    if not isinstance(A, np.ndarray):
        raise TypeError(f"Weighting matrix {name} must be a numpy array.")
    if A.ndim != 2:
        raise ValueError(f"Weighting matrix {name} must be a 2-dimensional matrix, got ndim={A.ndim}.")
    if A.shape[0] != A.shape[1]:
        raise ValueError(f"Weighting matrix {name} is not square.")
    if not np.allclose(A, A.T, atol=tol):
        raise warnings.warn(f"Weighting matrix {name} is not symmetric.")
    else:
        # check whether A is diagonal
        if np.all(np.abs(A - np.diag(np.diag(A))) < tol):
            if np.any(np.diag(A) < 0):
                raise warnings.warn(f"Diagonal weighting matrix {name} is not positive semi-definite.")
        else:
            try:
                E = np.linalg.eigvalsh(A)
            except:
                raise warnings.warn(f"Eigenvalue decomposition of weighting matrix {name} failed, the matrix might not be positive definite.")
            if not np.all(E > tol):
                raise warnings.warn(f"Weighting matrix {name} is not positive definite.")
