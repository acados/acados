
setenv("MW_MINGW64_LOC", "C:\ProgramData\mingw64\mingw64");
% setenv("MW_MINGW64_LOC", "C:\tools\mingw64\bin");
% setenv("MW_MINGW64_LOC",getenv("MW_MINGW64_LOC") + ";C:\ProgramData\mingw64\mingw64\bin");
mex -v -setup