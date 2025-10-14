function sol = moveInsert(P, Q, sol)
% 遗传算子变异操作
% 如果父母是SOLUTION类，则后代会是一个SOLUTION类；如果父母是决策变量，则后代也是决策变量
if isa(sol,'SOLUTION')
    calObj = true;
    sol = sol.decs;
else
    calObj = false;
end
if P < Q
    sol = sol([1:P-1, Q, P:Q-1, Q+1:end]);
elseif P > Q
    sol = sol([1:Q-1, Q+1:P-1, Q, P:end]);
end
if calObj
    sol = SOLUTION(sol);
end
end