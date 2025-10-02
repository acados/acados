function ns = ns_from_idxs_rev(idxs_rev)
    if isempty(idxs_rev)
        ns = 0;
    else
        ns = max(idxs_rev) + 1;
        for i=0:ns-1
            if ~ismember(i, idxs_rev)
                error(['Detected ns = ', num2str(ns), ', but i = ', num2str(i), ' is not in idxs_rev = [', num2str(idxs_rev), '], the slack with index ', num2str(i), ' is thus not contained in the problem.']);
            end
        end
    end
end
