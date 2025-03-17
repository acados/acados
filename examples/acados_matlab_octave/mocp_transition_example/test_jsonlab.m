% jsonlab test

addpath(fullfile('/home/oj/acados', 'external', 'jsonlab'));
json_filename = 'test.json';

x = struct();
x.a = 1;

in = {x, x}
json_string = savejson('', in, 'ForceRootName', 0)%, 'SimplifyCell',0, 'SingletCell', 0);
% fid = fopen(json_filename, 'w');
% if fid == -1, error('Cannot create json file'); end
% fwrite(fid, json_string, 'char');
% fclose(fid);
disp(json_string)
out = loadjson(json_string)%, 'SimplifyCell',0);
disp(out)
