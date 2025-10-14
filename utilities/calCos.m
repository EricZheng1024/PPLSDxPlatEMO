function cos = calCos(V1, V2)
% 计算向量夹角
    cos = 2 * ones(size(V1, 1), 1);
    if any(size(V1, 2) ~= size(V2, 2))
        return
    end
    V1V2Multiply = sum(V1.*V2, 2);
    V1Len = sqrt(sum(V1.^2, 2));
    V2Len = sqrt(sum(V2.^2, 2));
    index = V1Len ~= 0 & V2Len ~= 0;
    cos(index) = V1V2Multiply(index) ./ V1Len(index) ./ V2Len(index);
end


% function cos = calCos(V1, V2)
% % 计算向量夹角
%     if length(V1) ~= length(V2)
%         cos = 2;
%         return
%     end
%     V1V2Multiply = sum(V1.*V2);
%     V1Len = sqrt(sum(V1.^2));
%     V2Len = sqrt(sum(V2.^2));
%     if V1Len == 0 || V2Len == 0
%         cos = 2;
%     else
%         cos = V1V2Multiply / V1Len / V2Len;
%     end
% end