function m_mat = py2m(py_array, np)
%	dim = py_array.shape;
	ndim = int64(np.ndim(py_array));
	shape = np.shape(py_array);
%	py_item = py_array.item();
	if (ndim==1)
		m = int64(shape{1});
		m_mat = zeros(m, 1);
		for ii=1:m
			m_mat(ii) = py_array.item(int32(ii-1));
%			m_mat(ii) = item(int32(ii-1));
		end
	else % ndim==2
		m = int64(shape{1});
		n = int64(shape{2});
		m_mat = zeros(m, n);
		for ii=1:m
			for jj=1:n
				m_mat(ii, jj) = py_array.item(int32(ii-1), int32(jj-1));
%				m_mat(ii, jj) = item(int32(ii-1), int32(jj-1));
			end
		end
	end
end
