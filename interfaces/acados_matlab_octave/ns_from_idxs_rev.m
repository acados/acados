function n = ns_from_idxs_rev(idxs_rev)
    if isempty(idxs_rev)
        n = 0;
    else
        n = max(idxs_rev) + 1;
    end
end
