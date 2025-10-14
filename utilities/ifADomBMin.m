function res = ifADomBMin(obj1, obj2)
% 判断A是否支配B
if length(obj1) ~= length(obj2)
    res = false;
    return
end
res = all(obj1 <= obj2) && any(obj1 < obj2);
end