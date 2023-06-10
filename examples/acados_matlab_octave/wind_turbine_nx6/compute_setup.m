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

%% initial problem setup

nsim = 1000;

horLength = 40;
Ts = 0.2;
Tp = Ts*horLength;
N = ceil(Tp/Ts);

%% sim setup

% load wind field data
load ssOptimalTrackingRef.mat
ppOmegaSS = pchip(vwind,omegaSS);
ppColPitchSS = pchip(vwind,colPitchSS);

% Maximum wind speed for reference interpolation
VwindMaxRef = max(vwind);

% load FAST simulation data
load Sim_FAST__11mps_001_BC_FullSim.mat
% add inputs as additional 2 states
statesFAST(:,21:22) = uCtrlFAST(:,1:2);

% skip initial simulation data
tBeg = 40;
ind = find(tFAST-tBeg>0,1);

statesFAST  = statesFAST(ind:end,:);
uCtrlFAST   = uCtrlFAST(ind:end,:);
windFAST    = windFAST(ind:end);
tFAST       = tFAST(ind:end)-tFAST(ind);

XVgl        = statesFAST;
timeVgl     = tFAST;

tsim        = [0:Ts:floor(tFAST(end))];              % Simulation horizon
VwindSim    = interp1(tFAST,windFAST,tsim);      % Disturbance input signal

% Initial state of MPC
uCtrlFAST(1,:) = max([uCtrlFAST(1,:);zeros(size(uCtrlFAST(1,:)))]);
%                   GEN_agvelSt     DT_agvelTorsSt      GEN_agSt          DT_agTorsSt           BLD_agPtchActSt     GEN_trqActSt   BLD_agPtchDesSt      GEN_trqDesSt     
X0 = [statesFAST(1,[9               10                  19                20])                  uCtrlFAST(1,1)      uCtrlFAST(1,2) uCtrlFAST(1,1)       uCtrlFAST(1,2)]; 
X0([2 4]) = 0;
U0 = [0 0];   
X0(1) = X0(1);

%% references

%stage cost
y = [];
for i=1:nsim+N+1
    currVwind = VwindSim(i);
    currVwind = min(max(currVwind,0.1),VwindMaxRef);
    y = cat(1,y,[ppval(ppOmegaSS,currVwind) ppval(ppColPitchSS,currVwind) 0 0]);
end;



% setup references
x0_ref = X0';
u0_ref = U0';
wind0_ref = VwindSim(1:nsim);
y_ref = y(1:nsim,:)';

%size(x0_ref)
%size(u0_ref)
%size(wind0_ref)
%size(y_ref)


%% save to mat file
ny = 4;

tmp_windN_ref = [Ts*[0:nsim-1]; zeros(N+1,nsim)];
for ii=1:nsim
%	tmp_windN_ref(1+[1:N+1],ii) = VwindSim(ii)*ones(N+1,1);
	tmp_windN_ref(1+[1:N+1],ii) = VwindSim(ii+[0:N]);
end
save('windN_ref', 'tmp_windN_ref');

tmp_wind0_ref = [Ts*[0:nsim-1]; VwindSim(1:nsim)];
save('wind0_ref', 'tmp_wind0_ref');

%tmp_y_ref = [Ts*[1:nsim]; y_ref];
%save('y_ref', 'tmp_y_ref');
yt = y';
tmp_y_ref = [Ts*[0:nsim-1]; zeros(ny*N,nsim)];
for ii=1:nsim
	for jj=1:N
%		tmp_y_ref(1+(jj-1)*ny+[1:ny],ii) = yt(:,ii);
		tmp_y_ref(1+(jj-1)*ny+[1:ny],ii) = yt(:,ii+jj-1);
	end
end
save('y_ref', 'tmp_y_ref');

tmp_y_e_ref = [Ts*[0:nsim-1]; y_ref(1:2,:)];
save('y_e_ref', 'tmp_y_e_ref');



return;







%% print to file

% open file
myfile = fopen("setup.c", "w");

fprintf(myfile, "int nsim = %d;\n", nsim);

% x0_ref
fprintf(myfile, "\n");
fprintf(myfile, "double x0_ref[] = {\n");
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e, %1.15e\n", X0(1), X0(2), X0(3), X0(4), X0(5), X0(6), X0(7), X0(8));
fprintf(myfile, "};\n");

% u0
fprintf(myfile, "\n");
fprintf(myfile, "double u0_ref[] = {0, 0};\n");

% wind
fprintf(myfile, "\n");
fprintf(myfile, "double wind0_ref[] = {\n");
for ii=1:nsim
	fprintf(myfile, "%1.15e,\n", VwindSim(ii));
end
ii = nsim+1;
fprintf(myfile, "%1.15e\n", VwindSim(ii));
fprintf(myfile, "};\n");

% y
fprintf(myfile, "\n");
fprintf(myfile, "double y_ref[] = {\n");
for ii=1:nsim
	fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e,\n", y(ii,1), y(ii,2), y(ii,3), y(ii,4));
end
ii = nsim+1;
fprintf(myfile, "%1.15e, %1.15e, %1.15e, %1.15e\n", y(ii,1), y(ii,2), y(ii,3), y(ii,4));
fprintf(myfile, "};\n");

% close file
fclose(myfile);


