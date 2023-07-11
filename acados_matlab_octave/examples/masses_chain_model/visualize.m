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

% script to visualize the chain of masses

%drawnow update
drawnow('expose')

figure(1), set(gcf, 'Color','white');
clf

subplot(3,3,[1:6]);
tol = 0.00;
p = patch([-0.2, 1.2, 1.2, -0.2], [wall-tol, wall-tol, wall-tol, wall-tol], [-4, -4, 1, 1], 'g');
hold on;



tmp_pos = reshape(cur_pos, 3, length(cur_pos)/3);
cur_pos = zeros(3,1);
for ii=1:nfm
	cur_pos = [cur_pos, tmp_pos(:,1+2*(ii-1))];
end

plot3(cur_pos(1,:), cur_pos(2,:), cur_pos(3,:), '-ob', ...
'MarkerSize', 7.5, 'MarkerFaceColor', 'b', 'linewidth', 0.2);

view([-135 45*3/4]);

xlim([-0.2 1.2]);
ylim([-0.2 1.2]);
zlim([-4 1]);

grid on;

set(gca, 'Box', 'on');

xlabel( 'x [m]' );
ylabel( 'y [m]' );
zlabel( 'z [m]' );
