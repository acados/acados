%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function saveIntVector(v, vstr, directory)

fid = fopen([directory filesep vstr '.txt'], 'wt');

for i = 1:length(v)
    fprintf(fid,'%d\n',v(i));
end

fclose(fid);

end