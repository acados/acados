function ipiv = idx_perm_to_ipiv(idx_perm)

n = length(idx_perm);

vec = 1:n;

for ii=1:n
	idx0 = idx_perm(ii);
	for jj=ii:n
		if vec(jj)==idx0
			idx1 = jj;
			break
		end
	end
	tmp = vec(ii);
	vec(ii) =  vec(idx1);
	vec(idx1) = tmp;
	ipiv(ii) = idx1;
end
		
ipiv = ipiv-1; % C 0-based indexing
