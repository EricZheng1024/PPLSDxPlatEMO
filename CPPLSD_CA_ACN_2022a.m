function [A0_new, norm_scope, log_rho, n_restart] = CPPLSD_CA_ACN_2022a(problem, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same)
% ֧��2022a�İ汾����spmd��ص��滻Ϊlab�������޸ĵط�������"2022a"
% ����Լ����PPLS/D-C&A�㷨
% �ֻ�����״̬����̽��3*L_ACN���������󣬽���״̬ת��
% ����ʱ���Ծɽ�����Ŷ�������ģ�������PPLS/D-C&A��
% spmd
%   ͬһ�����ֵı�����spmd֮�󣬻���composite���ͣ�������ÿ��Ԫ�ؼ���Ӧworker�ĸ����ֱ���������ֵ������logs��A0_new������spmd���ʼ��
%   worker�����������Lһ��
%   ����ͬһ���ֱ����Ḵ��һ�ݽ���worker�����ʹ��problem_worker��ͳ��ÿ��worker�����۴�����parfor��ʵҲ���ƣ�
% δ��ʵ��C����Ϊ��
%   ׼��1*�������Ƿ��ڽ��ڣ���ʹ��WS�õ��Ľ������Ʈ������
%   ��ͬ�������rhoֵ��һ����ͬ����Ҫ�������ͨѶ
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

warning('�ѹ�ʱ�����������type_restart��')

n = length(A0(1).dec);  % ������������
m = length(A0(1).obj);  % Ŀ�����
cheRefPoint = max(A0.objs,[],1);  % �ο���
[~, L] = UniformPoint(L, m);  % �������������

% C&A���
% C
W = UniformPoint(L, m);  % ����Ȩ������
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
    l = labindex;

    problem_worker = problem;  % ��Ȼʵ����spmd�Ḵ��һ�ݣ�ֱ��PROBLEM.Current(problem)�ǿ��еģ����Է���spmd֮���problem.FE��1������������д������
    PROBLEM.Current(problem_worker);  % Ϊÿ��workers��ʼ������
    inBoundCount = 0;  % �����������ڸ���
    myWV = W(l,:);  % ��ǰ������Ȩ������
    myWV_ = 1./myWV ./ vecnorm(1./myWV,2,2);
    otherWV = W([1:l-1, l+1:end],:);  % ����������Ȩ������
    
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
            if judgeInBoundM2MMin(normalize_2(A0(i).obj, norm_scope), cheRefPoint, myWV, otherWV)
                inBoundCount = inBoundCount + 1;
                Al0 = [Al0, A0(i)];
            end
        end
        if inBoundCount == 0  % û�н����������ڣ��ҵ�������ĵ�
            % refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - normalize_2(A0.objs, norm_scope);
            [~, I] = sort(calCos(refPointMinusFit, repmat(myWV,length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end

    Al = Al0;  % ������û�н������������Ľ⼯
    explored = false(1, length(Al0));  % �Ƿ��ѽ�������������־
    timeUpFlag = false;  % ʱ���ֹ��־
    stepIndex = 0;  % 0�������ý�ֹͣ��1��recheck
    n_restart = 0;  % ͳ����������
    % start_time = tic;  % ��ʼʱ�䲻������ʼ��
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
                % weiFitness = min((cheRefPoint-Al(i_itSol).obj)./myWV);  % WeiFit_Chebyshev_Min
                % weiFitness = rho*sum((cheRefPoint-Al(i_itSol).obj).*myWV) + (1-rho)*min((cheRefPoint-Al(i_itSol).obj).*myWV_);  % CN
                weiFitness = rho*sum((cheRefPoint-normalize_2(Al(i_itSol).obj,norm_scope)).*myWV) + (1-rho)*min((cheRefPoint-normalize_2(Al(i_itSol).obj,norm_scope)).*myWV_);  % CN
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
                    % Nei = cell(n,n);
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
                        for k = 1 : n
                            % ��ʱ�ж�
                            % if toc(start_time) > maxRuntime
                            if cputime - start_time > maxRuntime
                                timeUpFlag = true;
                                break
                            end

                            % �����½�s'
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

                            % ACN���
                            count_c1_FE = count_c1_FE + 1;

                            if breakNeiExpFlag
                                break
                            end
                        end

                        % 2-bit flip
                        % for k1 = 1 : n - 1  % ������������
                        %     for k2 = k1 + 1 : n
                        %         % ��ʱ�ж�
                        %         % if toc(start_time) > maxRuntime
                        %         if cputime - start_time > maxRuntime
                        %             timeUpFlag = true;
                        %             break
                        %         end
                        % 
                        %         % �����½�s'
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
                        %         % ACN���
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
                        for k1 = 1 : n - 1  % ������������
                            for k2 = k1 + 1 : n
                                % ��ʱ�ж�
                                % if toc(start_time) > maxRuntime
                                if cputime - start_time > maxRuntime
                                    timeUpFlag = true;
                                    break
                                end
                                
                                % �����½�s'
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

            % C&A��� A
            % if ~isempty(tp) && toc(start_time) > maxRuntime * tp(1)
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
                            [tf,sourceOut] = labProbe('any');
                            if tf
                                E = [E labReceive(sourceOut)];
                                is_recv(sourceOut) = true;
                            end
                        end
                        E = E(NDSort(E.objs,1)==1);  % ������Լ���Ƿ�����ļ�飬���Ҫ��tpʱAl0�нⶼ����Լ��
                        norm_scope = [min(E.objs, [], 1); max(E.objs, [], 1)];
                        labSend(norm_scope, 2:L);
                    otherwise
                        labSend(El, 1);
                        norm_scope = labReceive(1);  % ��������
                end
                % ����
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
            
            % �����ǰ��s��Ȼ��Al0�У��������Ϊ�Ѿ�̽����
            for i_itSol = 1 : length(Al0)
                if judgeSameSolFit(Al0(i_itSol), sol)  % �ж�Al0(i_itSol)�ǲ���sol�����ܽ��ж�Ŀ�������Ƿ���ȣ�����i_itSol��ͣ������ȵĵ�һ����������������ѭ����ʹ��judgeSameSolFit�ܹ�����Ӧstrict_same�����浵�ж���Ŀ�꺯��ֵ��ͬ�Ľ�ʱ��judgeSameSolFit�����˻�ΪjudgeSameSolFit2
                    explored(i_itSol) = true;
                    break
                end
            end
            
            %����Al
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

            Al = Operator_GA4CO([Al0 Al0]);  % �Ŷ�����Ϊ���ظ��⣬Ĭ�ϲ�������
            explored = [true(1,length(Al0)) false(1,length(Al))];
            Al0 = [Al0  Al];

            n_restart = n_restart + 1;
        end
    end
    
    str = ['����'  num2str(n_restart)  '��'];
    % log = ['������'  char(string(l))  '������'  str  '������ʱ��'  num2str(toc(start_time))  's'];
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
    tmp = problem_worker{i};  % 2022aҪ����
    problem.FE = problem.FE + tmp.FE - 1;
end
PROBLEM.Current(problem);  % ע�͵�Ҳ��

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
