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

import os
import sys
import subprocess
import shutil
import urllib.request
from pathlib import Path
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build acados. "
                "Please install cmake and try again."
            )

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        
        # Create lib, include, and bin directories in the package
        package_dir = os.path.join(extdir, "acados_template")
        lib_dir = os.path.join(package_dir, "lib")
        include_dir = os.path.join(package_dir, "include")
        bin_dir = os.path.join(package_dir, "bin")
        os.makedirs(lib_dir, exist_ok=True)
        os.makedirs(include_dir, exist_ok=True)
        os.makedirs(bin_dir, exist_ok=True)
        
        # Use a temporary install directory
        install_dir = os.path.join(self.build_temp, "install")
        
        cmake_args = [
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_INSTALL_PREFIX={install_dir}",
        ]

        build_args = ["--config", "Release"]

        # Get the number of CPUs for parallel build
        if hasattr(os, "cpu_count"):
            build_args += ["-j", str(os.cpu_count())]
        else:
            build_args += ["-j4"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # Configure
        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        
        # Build
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )
        
        # Install to temporary directory
        subprocess.check_call(
            ["cmake", "--install", "."], cwd=self.build_temp
        )
        
        # Find the actual install directory (CMake might create nested directories)
        # Look for the lib directory
        actual_install_dir = install_dir
        if not os.path.exists(os.path.join(actual_install_dir, "lib")):
            # Check if it's nested in build/temp.../install
            nested_install = os.path.join(self.build_temp, "build", os.path.basename(self.build_temp), "install")
            if os.path.exists(nested_install):
                actual_install_dir = nested_install
        
        # Copy lib and include to package directory
        install_lib = os.path.join(actual_install_dir, "lib")
        install_include = os.path.join(actual_install_dir, "include")
        
        if os.path.exists(install_lib):
            for item in os.listdir(install_lib):
                src = os.path.join(install_lib, item)
                dst = os.path.join(lib_dir, item)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)
                elif os.path.isdir(src):
                    if os.path.exists(dst):
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
        
        if os.path.exists(install_include):
            for item in os.listdir(install_include):
                src = os.path.join(install_include, item)
                dst = os.path.join(include_dir, item)
                if os.path.isfile(src):
                    shutil.copy2(src, dst)
                elif os.path.isdir(src):
                    if os.path.exists(dst):
                        shutil.rmtree(dst)
                    shutil.copytree(src, dst)
        
        # Download and install tera renderer
        self.download_tera_renderer(bin_dir)
    
    def download_tera_renderer(self, bin_dir):
        """Download tera renderer binary during pip install"""
        import platform
        import urllib.request
        
        TERA_VERSION = "0.2.0"
        binary_ext = ".exe" if os.name == 'nt' else ""
        
        # Platform mapping
        PLATFORM2TERA = {
            'linux': 'linux',
            'darwin': 'osx',
            'win32': 'windows'
        }
        
        # Get architecture
        arch = platform.machine().lower()
        if arch in ['x86_64', 'amd64']:
            arch = 'amd64'
        elif arch in ['arm64', 'aarch64']:
            arch = 'arm64'
        else:
            print(f"Warning: Unsupported architecture {arch} for tera renderer. Skipping download.")
            print("You can manually download it from https://github.com/acados/tera_renderer/releases")
            return
        
        platform_name = PLATFORM2TERA.get(sys.platform)
        if not platform_name:
            print(f"Warning: Unsupported platform {sys.platform} for tera renderer. Skipping download.")
            print("You can manually download it from https://github.com/acados/tera_renderer/releases")
            return
        
        # Construct download URL
        repo_url = "https://github.com/acados/tera_renderer/releases"
        url = f"{repo_url}/download/v{TERA_VERSION}/t_renderer-v{TERA_VERSION}-{platform_name}-{arch}{binary_ext}"
        
        tera_path = os.path.join(bin_dir, f"t_renderer{binary_ext}")
        
        # Download tera
        try:
            print(f"Downloading tera renderer from {url}")
            with urllib.request.urlopen(url) as response, open(tera_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
            print("Successfully downloaded t_renderer.")
            
            # Make executable (Unix-like systems)
            if os.name != 'nt':
                os.chmod(tera_path, 0o755)
                print("Successfully made t_renderer executable.")
        except Exception as e:
            print(f"Warning: Failed to download tera renderer: {e}")
            print("You can manually download it from https://github.com/acados/tera_renderer/releases")
            print(f"Place it in: {bin_dir}")



setup(
    ext_modules=[CMakeExtension("acados_c")],
    cmdclass={"build_ext": CMakeBuild},
)
