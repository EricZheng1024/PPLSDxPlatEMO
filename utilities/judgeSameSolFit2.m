function res = judgeSameSolFit2(sol1, sol2)
% 判断两个解是否相同
% 对于TSP，代码的方法是先判断目标向量；若不相同，则判断决策向量是否等价（首先对两条路径分别找到其相同的节点索引，然后从该索引开始，根据进一步的判断正向或反向扫描）
% 该函数是与原方法等价的，甚至更具一般性，逻辑如下：
%   目标向量不相同，代表决策向量不相同
%   目标向量相同，不一定代表决策向量相同
if length(sol1.dec) ~= length(sol2.dec) || length(sol1.obj) ~= length(sol2.obj)
    error('决策向量或目标向量维度不相同。')
end
res = all(sol1.obj == sol2.obj);
end