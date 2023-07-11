%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

% load data
load u_x0.mat

% open file
myfile = fopen("u_x0.c", "w");

% x0
fprintf(myfile, "double x0[] = {\n%1.15e, %1.15e, %1.15e\n};\n", x0(1), x0(2), x0(3));

fprintf(myfile, "\n");

% u
fprintf(myfile, "int nsim = %d;\n", size(u,1));
fprintf(myfile, "\n");
fprintf(myfile, "double u_sim[] = {\n");
for ii=1:size(u,1)-1
	fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e,\n", u(ii,1), u(ii,2), u(ii,3), u(ii,4));
end
ii = size(u,1);
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e\n", u(ii,1), u(ii,2), u(ii,3), u(ii,4));
fprintf(myfile, "};\n");

% close file
fclose(myfile);
