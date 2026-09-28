%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% Licensed under the 2-Clause BSD License.

%

function saveDoubleMatrix(M, Mstr, directory)

fid = fopen([directory filesep Mstr '.txt'], 'wt');

for ii = 1:size(M,1)
    for jj = 1:size(M,2)
        fprintf(fid,'%1.16e\t',M(ii, jj));
    end
    fprintf(fid,'\n');
end

fclose(fid);

end