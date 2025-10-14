function [A0_new, norm_scope, log_rho, n_restart, log_tree] = CPPLSD_CA2_ACN(problem, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart, type_c1)
% ����Լ����PPLS/D-C&A�㷨
% ����CPPLSD_CA_ACN 2023.12.8ʵ��C��
%   ǿ��ʹ��׼��1*
%   ��������Ϊ������ʹ��������ͬ��rhoֵ
%   InitializeBestRecord��UpdateBestRecordû�п���Լ����֮�󽫿��ǽ�����Լ���ж�ɾ����������������ٶ�
% ɾ��һЩע��
% legacy��
% �ֻ�����״̬����̽��3*L_ACN���������󣬽���״̬ת��
% ����ʱ���Ծɽ�����Ŷ�������ģ�������PPLS/D-C&A��
% spmd
%   ͬһ�����ֵı�����spmd֮�󣬻���composite���ͣ�������ÿ��Ԫ�ؼ���Ӧworker�ĸ����ֱ���������ֵ������logs��A0_new������spmd���ʼ��
%   worker�����������Lһ��
%   ����ͬһ���ֱ����Ḵ��һ�ݽ���worker�����ʹ��problem_worker��ͳ��ÿ��worker�����۴�����parfor��ʵҲ���ƣ�
%   
%
% problem               platemo��������ʵ��
% A0                    ��֧��⼯���л��������飬����֧�䣬��ô֧����п��ܲ��ᱻ�滻
% L                     ����������
% maxRuntime            �������������ʱ��
% strict_same           ��ͬ��Ķ��壬0��ʾĿ��ֵ��ͬ��Ϊ��ͬ�⣬1��ʾ���߱�����ͬ������ͬ��
%
% ע�⣺
%   1.A0������������ĳ�ʼ��Ⱥ��
%   2.�����С�����⡣
%   3.֧�����У�permutation��������2-opt��inversion mutation����֧�ֶ����ƣ�binary��������1λ��ת��one bit flip��
%   4.�޻��ƹ���ͼ��
% 
% 2023.11.30  1. ��������ѡ��Al�����ŵĽ���Ϊ�����⵱ǰLS�Ķ���Ӧѡ�������⺯��ֵ��õģ���ԭ����CV_bestû�и���
%             2. ���ж����ʼ��ʱ���߽�������Ľ����ʼ������������������⺯��ֵ�ģ��ڽ��ڵĽ����ױ�������½�֧�䣬����update�����У���������н������Ľⲻ��������£�
%                Ҫ���������ұȵ�ǰ��������⺯��ֵ���õĽ��Ƿǳ��ѵģ���˻��ǽ���ʼ�⻮������
% 2023.12.1   ׼��1��2�����¼�����Ŀ�꺯��ֵ�������ظ����㣬������Nei�����������������
% 2023.12.5   1. ��ʼʱ�䲻������ʼ��
%             2. �Ľ�Al�ĸ���
%             3. ɾ��recheck��recheck��Ϊ������first move�Ľ⣩������ǰֹͣʱ�����н���Ϊδ��̽������ֹworkerͨ���������������ɾ��ǰ�İ汾���legacy����Ҫ���������ѭ�����������³�ʼ��������Ľ⣨������
% 2023.12.6   ʹ��CPUʱ��
% 2023.12.8   1. һ�������ɽ���������������㽫binary��permutation���ظ�����ϲ���������������һЩ�����������������ٶȣ�    ������������ڳ�ʼ��ʱ���ɵġ�����ɵ�˰��󡣡����������� 
%             2. �����������myMV��
% 2023.12.12  �޸��ο���δ���£�CPPLSD_CA_ACN�Ľ��ֻ�ܵ���PPLS/D����
% 2023.12.14  1. �޸�ACN��tau������£�ԭΪtau����̽����1����ʱ��1����Ϊtau��̽����3����ʱ��1
%             2. ���˳���������
% 2023.12.15  �����ʷ���յĽ�
% 2023.12.16  ��A��ͨ�Ÿ�tagΪ2����ֹ���յ�C���͵�SOLUTION������norm_scopeΪSOLUTION��
% 2023.12.18  �޸�g_cn��֧�ֶ����ͬʱ���㣨��������sum��min������ά�ȣ�


n = length(A0(1).dec);  % ������������
m = length(A0(1).obj);  % Ŀ�����
cheRefPoint = max(A0.objs,[],1);  % �ο���
[~, L] = UniformPoint(L, m);  % �������������

% C&A���
% C
W = UniformPoint(L, m);  % ����Ȩ������
W_ = 1./W ./ vecnorm(1./W,2,2);
B_tmp = pdist2(W,W,'cosine');
[~,B_tmp] = sort(B_tmp,2);
B_tmp = B_tmp(:,(1:K)+1);  % +1��Ϊ�˲������Լ�
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
    l = spmdIndex;
    rest = setdiff(1:L, l);  % ����������Ȩ����������

    warning('off', 'parallel:lang:mpi:IncomingMessageDiscarded')

    problem_worker = problem;  % ��Ȼʵ����spmd�Ḵ��һ�ݣ�ֱ��PROBLEM.Current(problem)�ǿ��еģ����Է���spmd֮���problem.FE��1������������д������
    PROBLEM.Current(problem_worker);  % Ϊÿ��workers��ʼ������
    inBoundCount = 0;  % �����������ڸ���
    
    % ACN���
    rho_cur = rho_init;
    state_rho = 1;
    tau = 1;
    Zeta = zeros(L_ACN, 3);
    tp_2 = 0 : 0.001 : 1;
    log_rho = zeros(size(tp_2));
    
    % Al0 = A0;  % �������ݻ��⼯
    
    % �ҵ����ڵ�ǰ������Ľ�
    if length(A0) == 1
        Al0 = A0;
    else
        Al0 = [];
        for i = 1 : length(A0)
            % if judgeInBoundM2MMin(A0(i).obj, cheRefPoint, myWV, otherWV)
            if judgeInBoundM2MMin(normalize_2(A0(i).obj, norm_scope), cheRefPoint, W(l,:), W(rest,:))
                inBoundCount = inBoundCount + 1;
                Al0 = [Al0, A0(i)];
            end
        end
        if inBoundCount == 0  % û�н����������ڣ��ҵ�������ĵ�
            % refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - normalize_2(A0.objs, norm_scope);
            [~, I] = sort(calCos(refPointMinusFit, repmat(W(l,:),length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end

    Al = Al0;  % ������û�н������������Ľ⼯
    explored = false(1, length(Al0));  % �Ƿ��ѽ�������������־
    timeUpFlag = false;  % ʱ���ֹ��־
    stepIndex = 0;  % 0�������ý�ֹͣ��1��recheck
    n_restart = 0;  % ͳ����������
    log_tree = Al0;

    % InitializeBestRecord    C&A��� C
    flag_updated = false(1, length(B{l}));
    if length(Al0) == 1
        solB = repmat(Al0, 1, length(B{l}));
    else
        solB = repmat(SOLUTION(), 1, length(B{l}));
        for i = 1 : length(B{l})
            g = g_cn(Al0.objs, cheRefPoint, W(B{l}(i),:), W_(B{l}(i),:), rho_cur);
            [~, I] = max(g);
            solB(i) = Al0(I);
        end
    end

    start_time = cputime;  % ��ʼʱ�䲻������ʼ��
    
    %% proceed
    while ~timeUpFlag
        while ~isempty(Al) && ~timeUpFlag
            % ACN���
            rho = min(max(rho_cur + (state_rho-2)*drho, 0), 1);
            count_c1_FE = 0;  % ׼��1���ĵ����۴���

            % �ҵ�������������ۺϺ���ֵ�Ľ⣨ע������С�����⣬�Ҳο����MOEA/D�Ĳ�ͬ��
            bestWeiFitness = -inf;
            CV_best = inf;  % Լ��Υ��ֵ
            for i_itSol = 1 : length(Al)
                weiFitness = g_cn(normalize_2(Al(i_itSol).obj,norm_scope), cheRefPoint, W(l,:), W_(l,:), rho);
                CV_Al = sum(max(0,Al(i_itSol).con));
                if CV_Al < CV_best || (weiFitness > bestWeiFitness && CV_Al == CV_best)  % ��Լ�������ж�Լ�������ж��ʶ�ֵ
                    bestWeiFitness = weiFitness;
                    CV_best = CV_Al;
                    bestItSol = Al(i_itSol);
                end
            end
            sol = bestItSol;
            
            % ����׼��
            acceptFlag = false;  % �Ƿ���ܵ�ǰ��
            switch problem.encoding  % �����������������
                case 'binary'
                    Nei = cell(1,n);
                case 'permutation'
                    Nei = cell(n,n);
            end
            for acceptCriIndex = 0 : 1
                if acceptCriIndex == 1 && acceptFlag == true
                    break
                end
                
                breakNeiExpFlag = false;  % ��������������׼��1ʱ��
                switch problem.encoding
                    case 'binary'
                        % for k = 1 : n
                        for k = randperm(n)
                            % ��ʱ�ж�
                            if cputime - start_time > maxRuntime
                                timeUpFlag = true;
                                break
                            end

                            % �����½�s'
                            if isempty(Nei{k})
                                solPrime = oneBitFlip(k, sol);
                                Nei{k} = solPrime;
                            else
                                solPrime = Nei{k};
                            end

                            % C&A��� C
                            % ע����rho_cur
                            if acceptCriIndex == 0
                                [solB, tmp] = UpdateBestRecord(solPrime, solB, cheRefPoint, W(B{l},:), W_(B{l},:), rho_cur);
                                flag_updated = flag_updated | tmp;
                            end

                            % C&A��� C
                            [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                                update_CA2_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                acceptCriIndex,solPrime,cheRefPoint,norm_scope,W(l,:),W_(l,:),rho,W(rest,:),strict_same,stepIndex,type_c1);  % in

                            if acceptFlag
                                log_tree = [log_tree solPrime];
                            end

                            count_c1_FE = count_c1_FE + 1;  % ACN���

                            if breakNeiExpFlag
                                break
                            end
                        end

                    case 'permutation'
                        % for k1 = 1 : n - 1  % ��������
                        %     for k2 = k1 + 1 : n
                        for k1 = randperm(n-1)  % ��������
                            for k2 = k1 + randperm(n-k1)
                                % ��ʱ�ж�
                                if cputime - start_time > maxRuntime
                                    timeUpFlag = true;
                                    break
                                end
                                
                                % �����½�s'
                                if isempty(Nei{k1,k2})

                                    switch class(problem)
                                        case 'MOTSP'
                                            solPrime = move2Opt(k1, k2, sol);
                                        case 'mQAP'
                                            solPrime = moveSwap(k1, k2, sol);
                                        otherwise
                                            solPrime = moveInsert(k1, k2, sol);
                                    end
                                    Nei{k1,k2} = solPrime;
                                else
                                    solPrime = Nei{k1,k2};
                                end

                                % C&A��� C
                                % ע����rho_cur
                                if acceptCriIndex == 0
                                    [solB, tmp] = UpdateBestRecord(solPrime, solB, cheRefPoint, W(B{l},:), W_(B{l},:), rho_cur);
                                    flag_updated = flag_updated | tmp;
                                end

                                % C&A��� C
                                [inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag] = ...
                                    update_CA2_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                    acceptCriIndex,solPrime,cheRefPoint,norm_scope,W(l,:),W_(l,:),rho,W(rest,:),strict_same,stepIndex,type_c1);  % in

                                if acceptFlag
                                    log_tree = [log_tree solPrime];
                                end

                                % ACN���
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

                % ACN���
                if acceptCriIndex == 0
                    Zeta(tau,state_rho) = count_c1_FE;  % ����tau��state_rho����tau��state_rhoǰ����
                    if state_rho == 3  % ����state_rho����state_rhoǰ����
                        tau = tau + 1;
                    end
                    state_rho = mod(state_rho, 3) + 1;
                    if tau > L_ACN
                        [~, I] = min(sum(Zeta, 1) - 0.1*(1:3));
                        rho_cur = min(max(rho_cur + (I-2)*drho, 0), 1);
                        tau = 1;
                        Zeta = zeros(L_ACN, 3);
                    end
                end
            end

            % C&A���
            % C
            % SendToNeighbors    ��΢��һ����"if the previous non-blocking sending procedure to process ki has been completed then"
            for i = 1 : length(B{l})
                if flag_updated(i)
                    spmdSend(solB(i), B{l}(i), 1);  % �������tag����ֹ��A��ͻ
                end
            end
            flag_updated = false(1, length(B{l}));
            % ReceiveFromNeighbors
            [tf,sourceOut] = spmdProbe('any', 1);
            while tf
                % Recv
                solRecv = spmdReceive(sourceOut, 1);
                % Update
                for i = 1 : length(solRecv)  % update_CA2_ACN��֧��һ���Ը��¶���⣬�ᱨ���������̫��
                    solPrime = solRecv(i);
                    % Update itself
                    acceptFlag = false;  % �Ƿ���ܵ�ǰ��
                    for acceptCriIndex = 0 : 1
                        if acceptCriIndex == 1 && acceptFlag == true
                            break
                        end
                        [inBoundCount,acceptFlag,Al0,explored,~] = ...
                            update_CA2_ACN(inBoundCount,acceptFlag,Al0,explored,[], ...  % in and out
                            acceptCriIndex,solPrime,cheRefPoint,norm_scope,W(l,:),W_(l,:),rho,W(rest,:),strict_same,stepIndex,type_c1);  % in
                    end
                    % Update neighbors
                    solB = UpdateBestRecord(solPrime, solB, cheRefPoint, W(B{l},:), W_(B{l},:), rho_cur);
                end
                % Probe
                [tf,sourceOut] = spmdProbe('any', 1);
            end
            % A  Ҫ��C֮��
            if ~isempty(tp) && cputime - start_time > maxRuntime * tp(1)
                tp(1) = [];
                [~, index_El_lb] = min(Al0.objs, [], 1);
                [~, index_El_ub] = max(Al0.objs, [], 1);
                El = Al0([index_El_lb index_El_ub]);
                switch l
                    case 1  % ��1Ϊ���߳�
                        E = El;
                        is_recv = false(1, L);
                        is_recv(1) = true;
                        while ~all(is_recv)
                            [tf,sourceOut] = spmdProbe('any', 2);
                            if tf
                                E = [E spmdReceive(sourceOut, 2)];
                                is_recv(sourceOut) = true;
                            end
                        end
                        E = E(NDSort(E.objs,1)==1);  % ������Լ���Ƿ�����ļ�飬���Ҫ��tpʱAl0�нⶼ����Լ��
                        norm_scope = [min(E.objs, [], 1); max(E.objs, [], 1)];
                        spmdSend(norm_scope, 2:L, 2);
                    otherwise
                        spmdSend(El, 1, 2);
                        norm_scope = spmdReceive(1, 2);  % ��������
                end
                % ����
                cheRefPoint = norm_scope(2,:);
                cheRefPoint = normalize_2(cheRefPoint, norm_scope);
                inBoundCount = 0;
                explored(:) = false;
                for i = 1 : length(Al0)
                    if judgeInBoundM2MMin(normalize_2(Al0(i).obj, norm_scope), cheRefPoint, W(l,:), W(rest,:))
                        inBoundCount = inBoundCount + 1;
                    end
                end
            end


            % 
            if ~isempty(tp_2) && cputime - start_time > maxRuntime * tp_2(1)
                log_rho(length(tp_2)) = rho_cur;
                tp_2(1) = [];
            end
            
            % �����ǰ��s��Ȼ��Al0�У��������Ϊ�Ѿ�̽����
            for i_itSol = 1 : length(Al0)
                if judgeSameSolFit(Al0(i_itSol), sol)  % �ж�Al0(i_itSol)�ǲ���sol�����ܽ��ж�Ŀ�������Ƿ���ȣ�����i_itSol��ͣ������ȵĵ�һ����������������ѭ����ʹ��judgeSameSolFit�ܹ�����Ӧstrict_same�����浵�ж���Ŀ�꺯��ֵ��ͬ�Ľ�ʱ��judgeSameSolFit�����˻�ΪjudgeSameSolFit2
                    explored(i_itSol) = true;
                    break
                end
            end
            
            %����Al
            Al = Al0(~explored);
        end
        
        % recheck
        if timeUpFlag
            break
        end
        if stepIndex == 0
            switch type_restart
                case 1  % ԭʼ����
                    Al = Al0;
                    explored = false(1, length(Al0));
                case 2
                    Al = Operator_GA4CO([Al0 Al0]);  % �Ŷ�����Ϊ���ظ��⣬Ĭ�ϲ�������
                    explored = [true(1,length(Al0)) false(1,length(Al))];
                    Al0 = [Al0  Al];
                case 3
                    if n_restart > 0
                        Al0_elite = [Al0 Al0_last];
                        Al0_elite = Al0_elite(NDSort(Al0_elite.objs, 1) == 1);  % ������Ƿ�����Լ��
                    else
                        Al0_elite = Al0;
                    end
                    Al0 = Operator_GA4CO([Al0_elite Al0_elite]);  % �Ŷ�����Ϊ���ظ��⣬Ĭ�ϲ�������
                    Al = Al0;
                    explored = false(1,length(Al0));
                    Al0_last = Al0;
                case 4
                    % ʹ�ý��棬���������Ŷ�
                    if n_restart > 0
                        Al0_elite = [Al0 Al0_last];
                        Al0_elite = Al0_elite(NDSort(Al0_elite.objs, 1) == 1);  % ������Ƿ�����Լ��
                    else
                        Al0_elite = Al0;
                    end
                    Al0 = Operator_GA4CO(Al0_elite(randi(length(Al0_elite),1,2*length(Al0_elite))));  % ���Ŷ�����3�����ﲻһ��
                    Al = Al0;
                    explored = false(1,length(Al0));
                    Al0_last = Al0;
            end

            n_restart = n_restart + 1;
        end
    end
    
    str = ['����'  num2str(n_restart)  '��'];
    log = ['������'  char(string(l))  '������'  str  '��CPUʱ��'  num2str(cputime - start_time)  's'];
    log_rho = fliplr(log_rho);
end

%%
for i = 1 : length(log)
    disp(char(log{i}));
end
disp(['��������������ʱ��'  num2str(toc(start_time_total))  's']);  % worker��cputime�Ƕ�����
A0_new = [];
for i = 1 : length(Al0)
    A0_new = [A0_new Al0{i}];
end
A0_new = A0_new(NDSort(A0_new.objs, A0_new.cons, inf) == 1);  % ������з�֧������
norm_scope = norm_scope{1};

for i = 1 : length(problem_worker)
    problem.FE = problem.FE + problem_worker{i}.FE - 1;
end
PROBLEM.Current(problem);  % ע�͵�Ҳ��

% Saving Composites is not supported!
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

tmp = log_tree;
log_tree = cell(1,L);
for i = 1 : L
    log_tree{i} = tmp{i};
end

end
