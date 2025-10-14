function A0_new = CPPLSD_ACN(problem, A0, L, maxRuntime, strict_same, is_recheck, rho_init, drho)
% 带有约束的PPLS/D算法
% 轮换采用状态，在探索3*L_ACN个解的邻域后，进行状态转移
%
% problem               platemo的问题类实例
% A0                    非支配解集，行或者列数组，若非支配，那么支配解有可能不会被替换
% L                     子区域数量
% maxRuntime            子区域最大运行时间
% strict_same           相同解的定义，0表示目标值相同即为相同解，1表示决策变量相同才是相同解
%
% 注意：
%   1.A0是任意子区域的初始种群。
%   2.求解最小化问题。
%   3.支持排列（permutation），采用2-opt（inversion mutation）；支持二进制（binary），采用1位翻转（one bit flip）
%   4.无绘制过程图。
% 
% 2023.11.30  1. 修正错误选择Al最大序号的解作为子问题当前LS的对象（应选择子问题函数值最好的），原因是CV_best没有更新
%             2. 当有多个初始解时，边界子问题的界外初始解往往是有最好子问题函数值的（在界内的解容易被界外的新解支配，但在update函数中，如果界内有解则界外的解不会产生更新）
%                要产生界内且比当前最好子问题函数值更好的解是非常难的，因此还是将初始解划分区域
% 2023.12.1   准则1和2都重新计算了目标函数值，存在重复计算，引入了Nei储存已评估的邻域解

A0_new = [];
n = length(A0(1).dec);  % 决策向量长度
m = length(A0(1).obj);  % 目标个数
cheRefPoint = max(A0.objs,[],1);  % 参考点
[~, L] = UniformPoint(L, m);  % 修正子区域个数

% ACN相关
L_ACN = 10;
% drho = 0.1;

logs = strings(1, L);

start_time_total = tic;
parfor l = 1 : L
    W = UniformPoint(L, m);  % 生成权重向量；此处生成W，防止W成为广播变量
    start_time = tic;
    PROBLEM.Current(problem);  % 为每个workers初始化问题
    inBoundCount = 0;  % 解在子区域内个数
    myWV = W(l,:);  % 当前子区域权重向量
    myWV_ = 1./myWV ./ vecnorm(1./myWV,2,2);
    otherWV = W([1:l-1, l+1:end],:);  % 其他子区域权重向量
    
    % ACN相关
    rho_cur = rho_init;
    state_rho = 1;
    tau = 1;
    Zeta = zeros(L_ACN, 3);
    
    % Al0 = A0;  % 子区域演化解集
    
    % 找到属于当前子区域的解
    if length(A0) == 1
        Al0 = A0;
    else
        Al0 = [];
        for i = 1 : length(A0)
            if judgeInBoundM2MMin(A0(i).obj, cheRefPoint, myWV, otherWV)
                inBoundCount = inBoundCount + 1;
                Al0 = [Al0, A0(i)];
            end
        end
        if inBoundCount == 0  % 没有解在子区域内，找到离最近的点
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            [~, I] = sort(calCos(refPointMinusFit, repmat(myWV,length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end

    Al = Al0;  % 子区域没有进行邻域搜索的解集
    explored = false(1, length(Al0));  % 是否已进行邻域搜索标志
    timeUpFlag = false;  % 时间截止标志
    
    for stepIndex = 0 : is_recheck  % 0：遇到好解停止；1：recheck
        while ~isempty(Al) && ~timeUpFlag
            % ACN相关
            rho = min(max(rho_cur + (state_rho-2)*drho, 0), 1);
            count_c1_FE = 0;  % 准则1消耗的评价次数

            % 找到子区域具有最大聚合函数值的解（注意是最小化问题，且参考点和MOEA/D的不同）
            bestWeiFitness = -inf;
            CV_best = inf;  % 约束违反值
            for i_itSol = 1 : length(Al)
                % weiFitness = min((cheRefPoint-Al(i_itSol).obj)./myWV);  % WeiFit_Chebyshev_Min
                weiFitness = rho*sum((cheRefPoint-Al(i_itSol).obj).*myWV) + (1-rho)*min((cheRefPoint-Al(i_itSol).obj).*myWV_);  % CN
                CV_Al = sum(max(0,Al(i_itSol).con));
                if CV_Al < CV_best || (weiFitness > bestWeiFitness && CV_Al == CV_best)  % 带约束，先判断约束，后判断适度值
                    bestWeiFitness = weiFitness;
                    CV_best = CV_Al;
                    bestItSol = Al(i_itSol);
                end
            end
            sol = bestItSol;
            
            % 接收准则
            acceptFlag = false;  % 是否接受当前解
            switch problem.encoding  % 储存已评估的邻域解
                case 'binary'
                    Nei = cell(1,n);
                case 'permutation'
                    Nei = cell(n,n);
            end
            for acceptCriIndex = 0 : 1
                if acceptCriIndex == 1 && acceptFlag == true
                    break
                end
                
                breakNeiExpFlag = false;  % 跳出邻域搜索，准则1时用
                switch problem.encoding
                    case 'binary'
                        for k = 1 : n
                            % 超时判断
                            if toc(start_time) > maxRuntime
                                timeUpFlag = true;
                                break
                            end
                            
                            % 产生新解s'
                            % solPrime = oneBitFlip(k, sol);
                            if isempty(Nei{k})
                                solPrime = oneBitFlip(k, sol);
                                Nei{k} = solPrime;
                            else
                                solPrime = Nei{k};
                            end
                            
                            [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                                update_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                acceptCriIndex,solPrime,cheRefPoint,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in
                            
                            % ACN相关
                            count_c1_FE = count_c1_FE + 1;

                            if breakNeiExpFlag
                                break
                            end
                        end
                        
                    case 'permutation'
                        for k1 = 1 : n - 1  % 遍历整个邻域
                            for k2 = k1 + 1 : n
                                % 超时判断
                                if toc(start_time) > maxRuntime
                                    timeUpFlag = true;
                                    break
                                end
                                
                                % 产生新解s'
                                % solPrime = move2Opt(k1, k2, sol);
                                if isempty(Nei{k1,k2})
                                    solPrime = move2Opt(k1, k2, sol);
                                    Nei{k1,k2} = solPrime;
                                else
                                    solPrime = Nei{k1,k2};
                                end
                                
                                [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                                    update_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                    acceptCriIndex,solPrime,cheRefPoint,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in
                                
                                % ACN相关
                                count_c1_FE = count_c1_FE + 1;

                                if breakNeiExpFlag
                                    break
                                end
                            end
                            if breakNeiExpFlag
                                break
                            end
                        end
                end

                % ACN相关
                if acceptCriIndex == 0
                    Zeta(tau,state_rho) = count_c1_FE;
                    state_rho = mod(state_rho, 3) + 1;
                    tau = tau + 1;
                    if tau > L_ACN
                        [~, I] = min(sum(Zeta, 1) - 0.1*(1:3));
                        rho_cur = min(max(rho_cur + (I-2)*drho, 0), 1);
                        tau = 1;
                        Zeta = zeros(L_ACN, 3);
                    end
                end
            end
            
            % 如果当前解s仍然在Al0中，将它标记为已经探索的
            for i_itSol = 1 : length(Al0)
                if judgeSameSolFit(Al0(i_itSol), sol)  % 判断Al0(i_itSol)是不是sol；不能仅判断目标向量是否相等，否则i_itSol会停留在相等的第一个索引处，出现死循环；使用judgeSameSolFit能够自适应strict_same，当存档中都是目标函数值不同的解时，judgeSameSolFit近似退化为judgeSameSolFit2
                    explored(i_itSol) = true;
                    break
                end
            end
            
            %更新Al
            Al = [];
            for i_itSol = 1 : length(Al0)
                if explored(i_itSol) == false
                    Al = [Al, Al0(i_itSol)];
                end
            end
        end
        
        % recheck
        if timeUpFlag
            break
        end
        if stepIndex == 0
            Al = Al0;
            explored = false(1, length(Al0));
        end
    end
    
    if stepIndex == 0
        str = '没有进行recheck';
    else
        if timeUpFlag == 0
            str = '已进行完整recheck';
        else
            str = '已进行部分recheck';
        end
    end
    logs(l) = ['子区域'  char(string(l))  '结束，'  str  '，用时'  num2str(toc(start_time))  's'];
    A0_new = [A0_new Al0];
end
for i = 1 : length(logs)
    disp(char(logs(i)));
end
disp(['搜索结束，总用时'  num2str(toc(start_time_total))  's']);
A0_new = A0_new(NDSort(A0_new.objs, A0_new.cons, inf) == 1);  % 总体进行非支配排序
end
