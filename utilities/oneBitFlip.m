function sol = oneBitFlip(bit, sol)
% one bit flip, also can be k-bit flip
% 如果父母是SOLUTION类，则后代会是一个SOLUTION类；如果父母是决策变量，则后代也是决策变量
if isa(sol,'SOLUTION')
    calObj = true;
    sol = sol.decs;
else
    calObj = false;
end
sol(bit) = ~sol(bit);
if calObj
    sol = SOLUTION(sol);
end
end