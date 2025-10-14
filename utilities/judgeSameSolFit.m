function res = judgeSameSolFit(sol1, sol2)
% 判断两个解是否相同（严格意义）
if length(sol1.dec) ~= length(sol2.dec) || length(sol1.obj) ~= length(sol2.obj)
    error('决策向量或目标向量维度不相同。')
end
if any(sol1.obj ~= sol2.obj)  % 目标向量比较以节省计算资源
    res = false;
    return
end
res = all(sol1.dec == sol2.dec);  % 决策向量比较
end