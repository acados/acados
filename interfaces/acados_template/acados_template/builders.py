import os
import sys
from subprocess import call


class CMakeBuilder:
    """
    Class to work with the `CMake` build system.
    """
    def __init__(self):
        self._source_dir = None  # private source directory, this is set to code_export_dir
        self.build_dir = 'build'
        """The CMake build directory; default: :py:code:`build`."""
        self._build_dir = None  # private build directory, usually rendered to abspath(build_dir)
        self.generator = None
        """Defines the generator, options can be found via `cmake --help` under 'Generator'. Linux default 'Unix Makefiles'; default value: `None`."""
        self.build_targets = None
        """A comma-separated list of the build targets, if `None` then all targets will be build; type: List of strings; default: `None`."""
        self.options_on = None
        """List of strings as CMake options which are translated to '-D Opt[0] -D Opt[1] ...'; default: `None`."""

    # Generate the command string for handling the cmake command.
    def get_cmd1_cmake(self):
        defines_str = ''
        if self.options_on is not None:
            defines_arr = [f' -D {opt}=ON' for opt in self.options_on]
            defines_str = ' '.join(defines_arr)
        generator_str = ''
        if self.generator is not None:
            generator_str = f' -G "{self.generator}"'
        return f'cmake -D CMAKE_INSTALL_PREFIX="{self._source_dir}"{defines_str}{generator_str} -Wdev -S "{self._source_dir}" -B "{self._build_dir}"'

    # Generate the command string for handling the build.
    def get_cmd2_build(self):
        cmd = f'cmake --build "{self._build_dir}" -j$(nproc)'
        if self.build_targets is not None:
            cmd += f' -t {self.build_targets}'
        return cmd

    # Generate the command string for handling the install command.
    def get_cmd3_install(self, ):
        return f'cmake --install "{self._build_dir}"'

    def exec(self, code_export_dir):
        """
        Execute the compilation using `CMake` with the given settings.
        """
        self._source_dir = code_export_dir
        self._build_dir = os.path.abspath(self.build_dir)
        try:
            os.mkdir(self._build_dir)
        except FileExistsError as e:
            pass

        try:
            os.chdir(self._build_dir)
            cmd_str = self.get_cmd1_cmake()
            # print(f'call("{cmd_str})"')
            retcode = call(cmd_str, shell=True)
            if retcode != 0:
                raise RuntimeError(f'CMake command "{cmd_str}" was terminated by signal {retcode}')
            cmd_str = self.get_cmd2_build()
            # print(f'call("{cmd_str}")')
            retcode = call(cmd_str, shell=True)
            if retcode != 0:
                raise RuntimeError(f'Build command "{cmd_str}" was terminated by signal {retcode}')
            cmd_str = self.get_cmd3_install()
            # print(f'call("{cmd_str}")')
            retcode = call(cmd_str, shell=True)
            if retcode != 0:
                raise RuntimeError(f'Install command "{cmd_str}" was terminated by signal {retcode}')
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)
        except Exception as e:
            print("Execution failed:", e, file=sys.stderr)
            exit(1)
