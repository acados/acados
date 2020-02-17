#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
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

import os, sys
import urllib.request
import shutil
import numpy as np
from casadi import SX, MX, DM

ALLOWED_CASADI_VERSIONS = ('3.5.1', '3.4.5', '3.4.0')
TERA_VERSION = "0.0.30"

def get_acados_path():
    ACADOS_PATH = os.environ.get('ACADOS_SOURCE_DIR')
    if not ACADOS_PATH:
        acados_template_path = os.path.dirname(os.path.abspath(__file__))
        acados_path = os.path.join(acados_template_path, '../../../')
        ACADOS_PATH = os.path.realpath(acados_path)
    return ACADOS_PATH

def get_tera_exec_path():
    ACADOS_PATH = get_acados_path()
    return os.path.join(ACADOS_PATH, 'bin/t_renderer')

platform2tera = {
    "linux": "linux",
    "darwin": "osx",
    "win32": "window.exe"
}


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
        else:
            return False
    elif x == None or x == []:
        return False
    else:
        raise Exception("is_column expects one of the following types: np.ndarray, casadi.MX, casadi.SX."
                        + " Got: " + str(type(x)))

def is_empty(x):
    if isinstance(x, (MX, SX, DM)):
        return x.is_empty()
    elif x == None or x == []:
        return True
    else:
        raise Exception("is_empty expects one of the following types: casadi.MX, casadi.SX, None, empty list."
                        + " Got: " + str(type(x)))

def casadi_length(x):
    if isinstance(x, (MX, SX, DM)):
        return int(np.prod(x.shape))
    else:
        raise Exception("casadi_length expects one of the following types: casadi.MX, casadi.SX."
                        + " Got: " + str(type(x)))


def get_tera():
    tera_path = get_tera_exec_path()
    acados_path = get_acados_path()

    if os.path.exists(tera_path) and os.access(tera_path, os.X_OK):
        return tera_path

    repo_url = "https://github.com/acados/tera_renderer/releases"
    url = "{}/download/v{}/t_renderer-v{}-{}".format(
        repo_url, TERA_VERSION, TERA_VERSION, platform2tera[sys.platform])

    manual_install = 'For manual installation follow these instructions:\n'
    manual_install += '1 Download binaries from {}\n'.format(url)
    manual_install += '2 Copy them in {}/bin\n'.format(acados_path)
    manual_install += '3 Strip the version and platform from the binaries: '
    manual_install += 'as t_renderer-v0.0.30-X -> t_renderer)\n'
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

    if input(msg) == 'y':
        print("Dowloading {}".format(url))
        with urllib.request.urlopen(url) as response, open(tera_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        print("Successfully downloaded t_renderer.")
        os.chmod(tera_path, 0o755)
        return tera_path

    msg_cancel = "\nYou cancelled automatic download.\n\n"
    msg_cancel += manual_install
    msg_cancel += "Once installed re-run your script.\n\n"
    print(msg_cancel)

    sys.exit(1)


def render_template(in_file, out_file, template_dir, json_path):
    cwd = os.getcwd()
    if not os.path.exists(template_dir):
        os.mkdir(template_dir)
    os.chdir(template_dir)

    tera_path = get_tera()

    # setting up loader and environment
    acados_path = os.path.dirname(os.path.abspath(__file__))

    template_glob = acados_path + '/c_templates_tera/*'
    acados_template_path = acados_path + '/c_templates_tera'

    # call tera as system cmd
    os_cmd = "{tera_path} '{template_glob}' '{in_file}' '{json_path}' '{out_file}'".format(
        tera_path=tera_path,
        template_glob=template_glob,
        json_path=json_path,
        in_file=in_file,
        out_file=out_file
    )
    status = os.system(os_cmd)
    if (status != 0):
        raise Exception('Rendering of {} failed! Exiting.\n'.format(in_file))

    os.chdir(cwd)


## Conversion functions
def np_array_to_list(np_array):
    if isinstance(np_array, (np.ndarray)):
        return np_array.tolist()
    elif isinstance(np_array, (SX)):
        return DM(np_array).full()
    elif isinstance(np_array, (DM)):
        return np_array.full()
    else:
        raise(Exception(
            "Cannot convert to list type {}".format(type(np_array))
        ))


def dict2json(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json(v)

        # v_type = str(type(v).__name__)
        # out_key = '__' + v_type + '__' + k.split('__', 1)[-1]
        out_key = k.split('__', 1)[-1]
        out[k.replace(k, out_key)] = v
    return out


def dict2json_layout(d):
    """ Convert dictionary containing the description of
    of the ocp_nlp to JSON format by stripping the
    property mangling and adding array dimension info.
    ALL items of type String will be converted
    to type ndarrray!

    Parameters
    ----------
    d : dict
        dictionary containing the description of
        the ocp_nlp.

    Returns
    ------
    out: dict
        postprocessed dictionary.
    """
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            v = dict2json_layout(v)

        v_type = str(type(v).__name__)
        if v_type == 'list':
            v_type = 'ndarray'

        # add array number of dimensions?
        # if v_type == 'ndarray':
        #     v_type = v_type + '_' + str(len(v.shape))
        out_key = k.split('__', 1)[-1]

        if isinstance(v, dict):
            out[k.replace(k, out_key)] = v
        else:
            out[k.replace(k, out_key)] = [v_type]

    return out


def json2dict(ocp_nlp, ocp_dims):
    # load JSON layout
    current_module = sys.modules[__name__]
    acados_path = os.path.dirname(current_module.__file__)
    with open(acados_path + '/acados_layout.json', 'r') as f:
        ocp_nlp_layout = json.load(f)

    out = json2dict_rec(ocp_nlp, ocp_dims, ocp_nlp_layout)
    return out


def json2dict_rec(ocp_nlp, ocp_dims, ocp_nlp_layout):
    """ convert ocp_nlp loaded JSON to dictionary. Mainly convert
    lists to arrays for easier handling.
    Parameters
    ---------
    ocp_nlp : dict
        dictionary loaded from JSON to be post-processed.

    ocp_dims : instance of AcadosOcpDims

    ocp_nlp_layout : dict
        acados ocp_nlp layout.

    Returns
    -------
    out : dict
        post-processed dictionary.
    """
    out = {}
    for k, v in ocp_nlp.items():
        if isinstance(v, dict):
            v = json2dict_rec(v, ocp_dims, ocp_nlp_layout[k])

        v_type__ = str(type(v).__name__)
        out_key = k.split('__', 1)[-1]
        v_type = out_key.split('__')[0]
        out_key = out_key.split('__', 1)[-1]
        if 'ndarray' in ocp_nlp_layout[k]:
            if isinstance(v, int) or isinstance(v, float):
                v = np.array([v])
        if (v_type == 'ndarray' or v_type__ == 'list') and (ocp_nlp_layout[k][0] != 'str'):
            dims_l = []
            dims_names = []
            dim_keys = ocp_nlp_layout[k][1]
            for item in dim_keys:
                dims_l.append(ocp_dims[item])
                dims_names.append(item)
            dims = tuple(dims_l)
            if v == []:
                # v = None
                try:
                    v = np.reshape(v, dims)
                except:
                    raise Exception('acados -- mismatching dimensions for field {0}. Provided data has dimensions {1}, while associated dimensions {2} are {3}'.format(out_key, [], dims_names, dims))
                # v = []
            else:
                v = np.array(v)
                v_dims = v.shape
                if dims !=v_dims:
                    raise Exception('acados -- mismatching dimensions for field {0}. Provided data has dimensions {1}, while associated dimensions {2} are {3}'.format(out_key, v_dims, dims_names, dims))
        out[k.replace(k, out_key)] = v
    return out


def J_to_idx(J):
    nrows = J.shape[0]
    idx = np.zeros((nrows, ))
    for i in range(nrows):
        this_idx = np.nonzero(J[i,:])[0]
        if len(this_idx) != 1:
            raise Exception('Invalid J matrix structure detected, must contain one nonzero element per row. Exiting.')
        if this_idx.size > 0 and J[i,this_idx[0]] != 1:
            raise Exception('J matrices can only contain 1s. Exiting.')
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
            raise Exception('J_to_idx_slack: Invalid J matrix. Exiting. Found more than one nonzero in row ' + str(i))
        if this_idx.size > 0 and J[i,this_idx[0]] != 1:
            raise Exception('J_to_idx_slack: J matrices can only contain 1s, got J(' + str(i) + ', ' + str(this_idx[0]) + ') = ' + str(J[i,this_idx[0]]) )
    if not i_idx == ncol:
            raise Exception('J_to_idx_slack: J must contain a 1 in every column!')
    return idx