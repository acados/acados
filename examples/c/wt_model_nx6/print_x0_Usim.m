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
load testSim.mat

% open file
myfile = fopen("u_x0.c", "w");

% x0
%fprintf(myfile, "double x0[] = {\n%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e\n};\n", statesFAST(1,1), statesFAST(1,2), statesFAST(1,3), statesFAST(1,4), statesFAST(1,5), statesFAST(1,6));

%fprintf(myfile, "\n");

% nsim
fprintf(myfile, "int nsim = %d;\n", size(Usim,1));
fprintf(myfile, "\n");

% u
fprintf(myfile, "double u_sim[] = {\n");
for ii=1:size(Usim,1)-1
	fprintf(myfile, "%1.15e, %1.15e,\n", Usim(ii,1), Usim(ii,2));
end
ii = size(Usim,1);
fprintf(myfile, "%1.15e, %1.15e\n", Usim(ii,1), Usim(ii,2));
fprintf(myfile, "};\n");
fprintf(myfile, "\n");

% p
fprintf(myfile, "double p_sim[] = {\n");
for ii=1:size(Usim,1)-1
	fprintf(myfile, "%1.15e,\n", Usim(ii,3));
end
ii = size(Usim,1);
fprintf(myfile, "%1.15e\n", Usim(ii,3));
fprintf(myfile, "};\n");
fprintf(myfile, "\n");

% x_ref
fprintf(myfile, "double x_ref[] = {\n");
for ii=1:size(Usim,1)-1
	fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e,\n", statesFAST(ii,1), statesFAST(ii,2), statesFAST(ii,3), statesFAST(ii,4), statesFAST(ii,5), statesFAST(ii,6));
end
ii = size(Usim,1);
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e\n", statesFAST(ii,1), statesFAST(ii,2), statesFAST(ii,3), statesFAST(ii,4), statesFAST(ii,5), statesFAST(ii,6));
fprintf(myfile, "};\n");
fprintf(myfile, "\n");

% close file
fclose(myfile);
