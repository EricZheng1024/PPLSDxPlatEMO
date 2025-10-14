function sol = move2Opt(P, Q, sol)
% 2-opt操作
% 如果父母是SOLUTION类，则后代会是一个SOLUTION类；如果父母是决策变量，则后代也是决策变量
if P >= Q  % 规定P在左边，Q在右边，否则报错
    error('P should be smaller than Q in 2-opt!')
end
if isa(sol,'SOLUTION')
    calObj = true;
    sol = sol.decs;
else
    calObj = false;
end
sol = sol([1:P-1, Q:-1:P, Q+1:end]);  % Q+1>n时，Q+1:end返回空
if calObj
    sol = SOLUTION(sol);
end
end