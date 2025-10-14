function [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
    update_CA2_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
    acceptCriIndex,solPrime,cheRefPoint,norm_scope,myWV,myWV_,rho,otherWV,strict_same,stepIndex,type_c1)  % in
% 相比于update_CA_ACN修改了准则1

% 当子区域有解，跳过不在子区域的解
% if inBoundCount > 0 && ~judgeInBoundM2MMin(normalize_2(solPrime.obj,norm_scope), cheRefPoint, myWV, otherWV)
% if acceptCriIndex == 1 && inBoundCount > 0 && ~judgeInBoundM2MMin(normalize_2(solPrime.obj,norm_scope), cheRefPoint, myWV, otherWV)  % 准则1*
%     return
% end
if inBoundCount > 0 && ~judgeInBoundM2MMin(normalize_2(solPrime.obj,norm_scope), cheRefPoint, myWV, otherWV)
    if type_c1  % 准则1
        return
    else  % 准则1*
        if acceptCriIndex == 1
            return
        end
    end
end


if acceptCriIndex == 0  % 论文中的接受准则1
    solPrimeWeiFit = rho*sum((cheRefPoint-normalize_2(solPrime.obj,norm_scope)).*myWV) + (1-rho)*min((cheRefPoint-normalize_2(solPrime.obj,norm_scope)).*myWV_);  % CN
    CV_Prime = sum(max(0,solPrime.con));
    betterWeiFitFlag = true;
    for i_itSol = 1 : length(Al0)  % 检查是否在存档Al0中存在聚合函数值更好的解
        solWeiFit = rho*sum((cheRefPoint-normalize_2(solPrime.obj,norm_scope)).*myWV) + (1-rho)*min((cheRefPoint-normalize_2(Al0(i_itSol).obj,norm_scope)).*myWV_);  % CN
        CV_Al0 = sum(max(0,Al0(i_itSol).con));
        if CV_Al0 < CV_Prime || (solWeiFit >= solPrimeWeiFit && CV_Al0 == CV_Prime)  % 带约束
            betterWeiFitFlag = false;
            break
        end
    end
    acceptFlag = betterWeiFitFlag;
else  % acceptCriIndex == 1  论文中准则2
    inArchFlag = false;  % 新解是否是A10中的解
    beDomdFlag = false;  % 新解是否被Al0中的解支配
    for i_itSol = 1 : length(Al0)
        if (strict_same && judgeSameSolFit(Al0(i_itSol), solPrime)) || ...
                (~strict_same && judgeSameSolFit2(Al0(i_itSol), solPrime))
            inArchFlag = true;
            break
        end
        if ifADomBMin(Al0(i_itSol).obj, solPrime.obj)
            beDomdFlag = true;
            break
        end
    end
    CV_Prime = sum(max(0,solPrime.con));
    CV_Al0 = sum(max(0,Al0(i_itSol).con));
    acceptFlag = ~inArchFlag && ...
        (CV_Al0 > CV_Prime || (~beDomdFlag && CV_Al0 == CV_Prime));  % 带约束
end

if acceptFlag
    % 删除存档中被新解s'支配的解
    res_index = true(1,length(Al0));
    for i_itSol = 1 : length(Al0)
        CV_Prime = sum(max(0,solPrime.con));
        CV_Al0 = sum(max(0,Al0(i_itSol).con));
        if CV_Prime < CV_Al0 || (ifADomBMin(solPrime.obj, Al0(i_itSol).obj) && CV_Prime == CV_Al0)  % 带约束
            % if judgeInBoundM2MMin(Al0(i_itSol).obj, cheRefPoint, myWV, otherWV)
            if judgeInBoundM2MMin(normalize_2(Al0(i_itSol).obj,norm_scope), cheRefPoint, myWV, otherWV)
                inBoundCount = inBoundCount - 1;
            end
            res_index(i_itSol) = false;
        end
    end
    Al0 = Al0(res_index);
    explored = explored(res_index);
    
    % 将新解s'加入
    Al0 = [Al0 solPrime];
    explored = [explored false];
    
    % 如果新解在子区域内，则更新子区域解个数
    % if judgeInBoundM2MMin(solPrime.obj, cheRefPoint, myWV, otherWV)
    if judgeInBoundM2MMin(normalize_2(solPrime.obj,norm_scope), cheRefPoint, myWV, otherWV)
        inBoundCount = inBoundCount + 1;
    end
    
    % 若在准则1就找到了好解，停止
    if acceptCriIndex == 0
        if stepIndex == 0
            breakNeiExpFlag = true;
            return
        end
    end
end
end