function inBound = judgeInBoundM2MMin(obj, refPoint, myWV, otherWV)
% 最小化问题，判断解是否在某子区域内
% ***判断临近的向量就行了，待优化
if length(obj) ~= length(refPoint)
    error('目标向量维度与参考点维度不一致。')
end

numOtherProcs = size(otherWV, 1);

inBound = true;
refPointMinusFit = refPoint - obj;  % 为正

myCos = calCos(refPointMinusFit, myWV);  % 向量角度范围为0~180度，cos随者角度增大而减小
if myCos == 2
    inBound = false;
    return
end

for otherIndex = 1 : numOtherProcs
    otherCos = calCos(refPointMinusFit, otherWV(otherIndex, :));  % 向量角度范围为0~180度
    if otherCos == 2
        inBound = false;
        return
    end
    if otherCos > myCos
        inBound = false;
        break
    end
end
end