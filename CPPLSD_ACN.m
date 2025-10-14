function A0_new = CPPLSD_ACN(problem, A0, L, maxRuntime, strict_same, is_recheck, rho_init, drho)
% ����Լ����PPLS/D�㷨
% �ֻ�����״̬����̽��3*L_ACN���������󣬽���״̬ת��
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

A0_new = [];
n = length(A0(1).dec);  % ������������
m = length(A0(1).obj);  % Ŀ�����
cheRefPoint = max(A0.objs,[],1);  % �ο���
[~, L] = UniformPoint(L, m);  % �������������

% ACN���
L_ACN = 10;
% drho = 0.1;

logs = strings(1, L);

start_time_total = tic;
parfor l = 1 : L
    W = UniformPoint(L, m);  % ����Ȩ���������˴�����W����ֹW��Ϊ�㲥����
    start_time = tic;
    PROBLEM.Current(problem);  % Ϊÿ��workers��ʼ������
    inBoundCount = 0;  % �����������ڸ���
    myWV = W(l,:);  % ��ǰ������Ȩ������
    myWV_ = 1./myWV ./ vecnorm(1./myWV,2,2);
    otherWV = W([1:l-1, l+1:end],:);  % ����������Ȩ������
    
    % ACN���
    rho_cur = rho_init;
    state_rho = 1;
    tau = 1;
    Zeta = zeros(L_ACN, 3);
    
    % Al0 = A0;  % �������ݻ��⼯
    
    % �ҵ����ڵ�ǰ������Ľ�
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
        if inBoundCount == 0  % û�н����������ڣ��ҵ�������ĵ�
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            [~, I] = sort(calCos(refPointMinusFit, repmat(myWV,length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end

    Al = Al0;  % ������û�н������������Ľ⼯
    explored = false(1, length(Al0));  % �Ƿ��ѽ�������������־
    timeUpFlag = false;  % ʱ���ֹ��־
    
    for stepIndex = 0 : is_recheck  % 0�������ý�ֹͣ��1��recheck
        while ~isempty(Al) && ~timeUpFlag
            % ACN���
            rho = min(max(rho_cur + (state_rho-2)*drho, 0), 1);
            count_c1_FE = 0;  % ׼��1���ĵ����۴���

            % �ҵ�������������ۺϺ���ֵ�Ľ⣨ע������С�����⣬�Ҳο����MOEA/D�Ĳ�ͬ��
            bestWeiFitness = -inf;
            CV_best = inf;  % Լ��Υ��ֵ
            for i_itSol = 1 : length(Al)
                % weiFitness = min((cheRefPoint-Al(i_itSol).obj)./myWV);  % WeiFit_Chebyshev_Min
                weiFitness = rho*sum((cheRefPoint-Al(i_itSol).obj).*myWV) + (1-rho)*min((cheRefPoint-Al(i_itSol).obj).*myWV_);  % CN
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
                        for k = 1 : n
                            % ��ʱ�ж�
                            if toc(start_time) > maxRuntime
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
                                update_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                acceptCriIndex,solPrime,cheRefPoint,myWV,myWV_,rho,otherWV,strict_same,stepIndex);  % in
                            
                            % ACN���
                            count_c1_FE = count_c1_FE + 1;

                            if breakNeiExpFlag
                                break
                            end
                        end
                        
                    case 'permutation'
                        for k1 = 1 : n - 1  % ������������
                            for k2 = k1 + 1 : n
                                % ��ʱ�ж�
                                if toc(start_time) > maxRuntime
                                    timeUpFlag = true;
                                    break
                                end
                                
                                % �����½�s'
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
            
            % �����ǰ��s��Ȼ��Al0�У��������Ϊ�Ѿ�̽����
            for i_itSol = 1 : length(Al0)
                if judgeSameSolFit(Al0(i_itSol), sol)  % �ж�Al0(i_itSol)�ǲ���sol�����ܽ��ж�Ŀ�������Ƿ���ȣ�����i_itSol��ͣ������ȵĵ�һ����������������ѭ����ʹ��judgeSameSolFit�ܹ�����Ӧstrict_same�����浵�ж���Ŀ�꺯��ֵ��ͬ�Ľ�ʱ��judgeSameSolFit�����˻�ΪjudgeSameSolFit2
                    explored(i_itSol) = true;
                    break
                end
            end
            
            %����Al
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
        str = 'û�н���recheck';
    else
        if timeUpFlag == 0
            str = '�ѽ�������recheck';
        else
            str = '�ѽ��в���recheck';
        end
    end
    logs(l) = ['������'  char(string(l))  '������'  str  '����ʱ'  num2str(toc(start_time))  's'];
    A0_new = [A0_new Al0];
end
for i = 1 : length(logs)
    disp(char(logs(i)));
end
disp(['��������������ʱ'  num2str(toc(start_time_total))  's']);
A0_new = A0_new(NDSort(A0_new.objs, A0_new.cons, inf) == 1);  % ������з�֧������
end
