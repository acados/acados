%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function saveDoubleVector(v, vstr, directory)

fid = fopen([directory filesep vstr '.txt'], 'wt');

v = v(:);
for i = 1:length(v)
    fprintf(fid,'%1.16e\n',v(i));
end

fclose(fid);

end