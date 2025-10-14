function objs_n = normalize_2(objs, norm_scope)
    if isempty(norm_scope)
        objs_n = objs;
    else
        objs_n = (objs - norm_scope(1,:)) ./ (norm_scope(2,:) - norm_scope(1,:));
    end
end