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

ALLOWED_CASADI_VERSIONS = ('3.5.1', '3.4.5', '3.4.0')
TERA_VERSION = "0.0.30"

ACADOS_PATH = os.environ.get('ACADOS_SOURCE_DIR')
if not ACADOS_PATH:
    acados_template_path = os.path.dirname(os.path.abspath(__file__))
    acados_path = os.path.join(acados_template_path, '../../../')
    ACADOS_PATH = os.path.realpath(acados_path)

TERA_EXEC_PATH = os.path.join(ACADOS_PATH, 'bin/t_renderer')

platform2tera = {
    "linux": "linux",
    "darwin": "osx",
    "win32": "window.exe"
}

def is_column(x):
    if x.ndim == 1:
        return True
    elif x.ndim == 2 and x.shape[1] == 1:
        return True
    else:
        return False


def get_tera():
    tera_path = TERA_EXEC_PATH

    if os.path.exists(tera_path) and os.access(tera_path, os.X_OK):
        return tera_path

    repo_url = "https://github.com/acados/tera_renderer/releases"
    url = "{}/download/v{}/t_renderer-v{}-{}".format(
        repo_url, TERA_VERSION, TERA_VERSION, platform2tera[sys.platform])

    manual_install = 'For manual installation follow these instructions:\n'
    manual_install += '1 Download binaries from {}\n'.format(url)
    manual_install += '2 Copy them in {}/bin\n'.format(ACADOS_PATH)
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
