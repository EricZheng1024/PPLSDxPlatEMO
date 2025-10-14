function [A0_new, norm_scope, log_rho, n_restart] = CPPLSD_CA_ACN_2022a(problem, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same)
% 支持2022a的版本，将spmd相关的替换为lab，其他修改地方请搜索"2022a"
% 带有约束的PPLS/D-C&A算法
% 轮换采用状态，在探索3*L_ACN个解的邻域后，进行状态转移
% 重启时，对旧解进行扰动（额外的，不属于PPLS/D-C&A）
% spmd
%   同一个名字的变量在spmd之后，会变成composite类型，变量中每个元素即对应worker的该名字变量的最终值，所以logs和A0_new无需在spmd外初始化
%   worker的数量必须和L一致
%   由于同一名字变量会复制一份进入worker，因此使用problem_worker来统计每个worker的评价次数（parfor其实也类似）
% 未能实现C，因为：
%   准则1*不考虑是否在界内，这使得WS得到的解很容易飘到界外
%   不同子区域的rho值不一定相同，需要对其进行通讯
%   
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
% 2023.12.5   1. 开始时间不包含初始化
%             2. 改进Al的更新
%             3. 删除recheck（recheck是为了搜完first move的解），当提前停止时将所有解设为未被探索，防止worker通信造成永久阻塞，删除前的版本请见legacy，主要改了最外层循环；引入重新初始化子区域的解（待定）
% 2023.12.6   使用CPU时间

warning('已过时。请引入参数type_restart。')

n = length(A0(1).dec);  % 决策向量长度
m = length(A0(1).obj);  % 目标个数
cheRefPoint = max(A0.objs,[],1);  % 参考点
[~, L] = UniformPoint(L, m);  % 修正子区域个数

% C&A相关
% C
W = UniformPoint(L, m);  % 生成权重向量
B_tmp = pdist2(W,W,'cosine');
[~,B_tmp] = sort(B_tmp,2);
B_tmp = B_tmp(:,(1:K)+1);  % +1是为了不包含自己
B = cell(L,1);  % mutual K-closest
for i = 1 : L
    for j = 1 : K
        if ismember(i, B_tmp(B_tmp(i,j),:))
            B{i} = [B{i} B_tmp(i,j)];
        end
    end
end
% A
norm_scope = [];
cheRefPoint = normalize_2(cheRefPoint, norm_scope);

poolobj = gcp('nocreate');
if isempty(poolobj)  ||  poolobj.NumWorkers ~= L
    delete(poolobj);
    parpool(L);
end

%%
start_time_total = tic;
spmd
    %% Config
    l = labindex;

    problem_worker = problem;  % 虽然实际上spmd会复制一份，直接PROBLEM.Current(problem)是可行的（可以发现spmd之后的problem.FE是1），但是这样写更清晰
    PROBLEM.Current(problem_worker);  % 为每个workers初始化问题
    inBoundCount = 0;  % 解在子区域内个数
    myWV = W(l,:);  % 当前子区域权重向量
    myWV_ = 1./myWV ./ vecnorm(1./myWV,2,2);
    otherWV = W([1:l-1, l+1:end],:);  % 其他子区域权重向量
    
    % ACN相关
    rho_cur = rho_init;
    state_rho = 1;
    tau = 1;
    Zeta = zeros(L_ACN, 3);
    tp_2 = 0 : 0.001 : 1;
    log_rho = zeros(size(tp_2));
    
    % Al0 = A0;  % 子区域演化解集
    
    % 找到属于当前子区域的解
    if length(A0) == 1
        Al0 = A0;
    else
        Al0 = [];
        for i = 1 : length(A0)
            % if judgeInBoundM2MMin(A0(i).obj, cheRefPoint, myWV, otherWV)
            if judgeInBoundM2MMin(normalize_2(A0(i).obj, norm_scope), cheRefPoint, myWV, otherWV)
                inBoundCount = inBoundCount + 1;
                Al0 = [Al0, A0(i)];
            end
        end
        if inBoundCount == 0  % 没有解在子区域内，找到离最近的点
            % refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - normalize_2(A0.objs, norm_scope);
            [~, I] = sort(calCos(refPointMinusFit, repmat(myWV,length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end

    Al = Al0;  % 子区域没有进行邻域搜索的解集
    explored = false(1, length(Al0));  % 是否已进行邻域搜索标志
    timeUpFlag = false;  % 时间截止标志
    stepIndex = 0;  % 0：遇到好解停止；1：recheck
    n_restart = 0;  % 统计重启次数
    % start_time = tic;  % 开始时间不包含初始化
    start_time = cputime;  % 开始时间不包含初始化
    
    
    %% proceed
    while ~timeUpFlag
        while ~isempty(Al) && ~timeUpFlag
            % ACN相关
            rho = min(max(rho_cur + (state_rho-2)*drho, 0), 1);
            count_c1_FE = 0;  % 准则1消耗的评价次数

            % 找到子区域具有最大聚合函数值的解（注意是最小化问题，且参考点和MOEA/D的不同）
            bestWeiFitness = -inf;
            CV_best = inf;  % 约束违反值
            for i_itSol = 1 : length(Al)
                % weiFitness = min((cheRefPoint-Al(i_itSol).obj)./myWV);  % WeiFit_Chebyshev_Min
                % weiFitness = rho*sum((cheRefPoint-Al(i_itSol).obj).*myWV) + (1-rho)*min((cheRefPoint-Al(i_itSol).obj).*myWV_);  % CN
                weiFitness = rho*sum((cheRefPoint-normalize_2(Al(i_itSol).obj,norm_scope)).*myWV) + (1-rho)*min((cheRefPoint-normalize_2(Al(i_itSol).obj,norm_scope)).*myWV_);  % CN
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
                    % Nei = cell(n,n);
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
                            % if toc(start_time) > maxRuntime
                            if cputime - start_time > maxRuntime
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
                                update_CA_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                acceptCriIndex,solPrime,cheRefPoint,norm_scope,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in

                            % ACN相关
                            count_c1_FE = count_c1_FE + 1;

                            if breakNeiExpFlag
                                break
                            end
                        end

                        % 2-bit flip
                        % for k1 = 1 : n - 1  % 遍历整个邻域
                        %     for k2 = k1 + 1 : n
                        %         % 超时判断
                        %         % if toc(start_time) > maxRuntime
                        %         if cputime - start_time > maxRuntime
                        %             timeUpFlag = true;
                        %             break
                        %         end
                        % 
                        %         % 产生新解s'
                        %         if isempty(Nei{k1,k2})
                        %             solPrime = oneBitFlip([k1,k2], sol);
                        %             Nei{k1,k2} = solPrime;
                        %         else
                        %             solPrime = Nei{k1,k2};
                        %         end
                        % 
                        %         [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                        %             update_CA_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                        %             acceptCriIndex,solPrime,cheRefPoint,norm_scope,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in
                        % 
                        %         % ACN相关
                        %         count_c1_FE = count_c1_FE + 1;
                        % 
                        %         if breakNeiExpFlag
                        %             break
                        %         end
                        %     end
                        %     if breakNeiExpFlag
                        %         break
                        %     end
                        % end
                        
                    case 'permutation'
                        for k1 = 1 : n - 1  % 遍历整个邻域
                            for k2 = k1 + 1 : n
                                % 超时判断
                                % if toc(start_time) > maxRuntime
                                if cputime - start_time > maxRuntime
                                    timeUpFlag = true;
                                    break
                                end
                                
                                % 产生新解s'
                                % solPrime = move2Opt(k1, k2, sol);
                                if isempty(Nei{k1,k2})

                                    switch class(problem)
                                        case 'MOTSP'
                                            solPrime = move2Opt(k1, k2, sol);
                                        case 'mQAP'
                                            solPrime = moveSwap(k1, k2, sol);
                                        otherwise
                                            error('Operator undefined.')
                                    end
                                    Nei{k1,k2} = solPrime;
                                else
                                    solPrime = Nei{k1,k2};
                                end
                                
                                [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                                    update_CA_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                    acceptCriIndex,solPrime,cheRefPoint,norm_scope,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in
                                
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

            % C&A相关 A
            % if ~isempty(tp) && toc(start_time) > maxRuntime * tp(1)
            if ~isempty(tp) && cputime - start_time > maxRuntime * tp(1)
                tp(1) = [];
                [~, index_El_lb] = min(Al0.objs, [], 1);
                [~, index_El_ub] = max(Al0.objs, [], 1);
                El = Al0([index_El_lb index_El_ub]);
                switch l
                    case 1  % 设1为主线程
                        E = El;
                        is_recv = false(1, L);
                        is_recv(1) = true;
                        while ~all(is_recv)
                            [tf,sourceOut] = labProbe('any');
                            if tf
                                E = [E labReceive(sourceOut)];
                                is_recv(sourceOut) = true;
                            end
                        end
                        E = E(NDSort(E.objs,1)==1);  % 不进行约束是否满足的检查，因此要求tp时Al0中解都满足约束
                        norm_scope = [min(E.objs, [], 1); max(E.objs, [], 1)];
                        labSend(norm_scope, 2:L);
                    otherwise
                        labSend(El, 1);
                        norm_scope = labReceive(1);  % 阻塞接收
                end
                % 重置
                cheRefPoint = normalize_2(cheRefPoint, norm_scope);
                inBoundCount = 0;
                explored(:) = false;
                for i = 1 : length(Al0)
                    if judgeInBoundM2MMin(normalize_2(Al0(i).obj, norm_scope), cheRefPoint, myWV, otherWV)
                        inBoundCount = inBoundCount + 1;
                    end
                end
            end

            % 
            % if ~isempty(tp_2) && toc(start_time) > maxRuntime * tp_2(1)
            if ~isempty(tp_2) && cputime - start_time > maxRuntime * tp_2(1)
                log_rho(length(tp_2)) = rho_cur;
                tp_2(1) = [];
            end
            
            % 如果当前解s仍然在Al0中，将它标记为已经探索的
            for i_itSol = 1 : length(Al0)
                if judgeSameSolFit(Al0(i_itSol), sol)  % 判断Al0(i_itSol)是不是sol；不能仅判断目标向量是否相等，否则i_itSol会停留在相等的第一个索引处，出现死循环；使用judgeSameSolFit能够自适应strict_same，当存档中都是目标函数值不同的解时，judgeSameSolFit近似退化为judgeSameSolFit2
                    explored(i_itSol) = true;
                    break
                end
            end
            
            %更新Al
            % Al = [];
            % for i_itSol = 1 : length(Al0)
            %     if explored(i_itSol) == false
            %         Al = [Al, Al0(i_itSol)];
            %     end
            % end
            Al = Al0(~explored);
        end
        
        % recheck
        if timeUpFlag
            break
        end
        if stepIndex == 0
            % Al = Al0;
            % explored = false(1, length(Al0));

            Al = Operator_GA4CO([Al0 Al0]);  % 扰动；因为是重复解，默认参数即可
            explored = [true(1,length(Al0)) false(1,length(Al))];
            Al0 = [Al0  Al];

            n_restart = n_restart + 1;
        end
    end
    
    str = ['重启'  num2str(n_restart)  '次'];
    % log = ['子区域'  char(string(l))  '结束，'  str  '，挂钟时间'  num2str(toc(start_time))  's'];
    log = ['子区域'  char(string(l))  '结束，'  str  '，CPU时间'  num2str(cputime - start_time)  's'];
    log_rho = fliplr(log_rho);
end

%%
for i = 1 : length(log)
    disp(char(log{i}));
end
disp(['搜索结束，挂钟时间'  num2str(toc(start_time_total))  's']);  % worker的cputime是独立的
A0_new = [];
for i = 1 : length(Al0)
    A0_new = [A0_new Al0{i}];
end
A0_new = A0_new(NDSort(A0_new.objs, A0_new.cons, inf) == 1);  % 总体进行非支配排序
norm_scope = norm_scope{1};

for i = 1 : length(problem_worker)
    tmp = problem_worker{i};  % 2022a要这样
    problem.FE = problem.FE + tmp.FE - 1;
end
PROBLEM.Current(problem);  % 注释掉也行

tmp = log_rho;
log_rho = cell(1,L);
for i = 1 : L
    log_rho{i} = tmp{i};
end

tmp = n_restart;
n_restart = cell(1,L);
for i = 1 : L
    n_restart{i} = tmp{i};
end

end
