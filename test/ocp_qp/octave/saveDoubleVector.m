function saveDoubleVector(v, vstr, directory)

fid = fopen([directory filesep vstr '.txt'], 'wt');

v = v(:);
for i = 1:length(v)
    fprintf(fid,'%1.16e\n',v(i));
end

fclose(fid);

end