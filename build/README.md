In order to compile for dSPACE mex binaries are needed for the HOST machine as well. In order
to obtain them use CMake from Windows's Powershell (other soultions might be possible, 
this is what we used in practice):

'cmake -G "Visual Studio 15 2017 Win64" -DA
CADOS_WITH_OSQP=OFF -DACADOS_WITH_QPDUNES=OFF -DACADOS_WITH_HPMPC=OFF -DCMAKE_IN
STALL_PREFIX="..\acados_install_dir" ..'

and then 

'devenv /build Release acados.sln'