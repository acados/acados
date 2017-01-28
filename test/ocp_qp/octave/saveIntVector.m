function saveIntVector(v, vstr, directory)

fid = fopen([directory filesep vstr '.txt'], 'wt');

for i = 1:length(v)
    fprintf(fid,'%d\n',v(i));
end

fclose(fid);

end