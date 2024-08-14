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


function render_file( json_fullfile, template_dir, template_file, out_file, ...
                      t_renderer_location )

    os_cmd = [t_renderer_location, ' "',...
        template_dir, '"', ' ', '"', template_file, '"', ' ', '"',...
        json_fullfile, '"', ' ', '"', out_file, '"'];

    [ status, result ] = system(os_cmd);
    if status
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result:\n%s\n\n',...
            template_file, os_cmd, status, result);
    end
    % NOTE: this should return status != 0, maybe fix in tera renderer?
    if ~isempty(strfind( result, 'Error' )) % contains not implemented in Octave
        cd ..
        error('rendering %s failed.\n command: %s\n returned status %d, got result: %s',...
            template_file, os_cmd, status, result);
    end
end

