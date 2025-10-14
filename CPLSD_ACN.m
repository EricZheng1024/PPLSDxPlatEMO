function A0_new = CPLSD_ACN(problem, A0, L, maxRuntime, process_visualization, strict_same, is_recheck, rho_init, drho)
% ����Լ����PLS/D�㷨��û�в��е�PPLS/D��
% ������ACN���ԣ�Ҫ���ٶ�λ�޸�������"ACN���"
% �ֻ�����״̬����̽��3*L_ACN���������󣬽���״̬ת��
%
% problem               platemo��������ʵ��
% A0                    ��֧��⼯���л��������飬����֧�䣬��ô֧����п��ܲ��ᱻ�滻
% L                     ����������
% maxRuntime            �������������ʱ��
% process_visualization �Ƿ���й��̿��ӻ�
% strict_same           ��ͬ��Ķ��壬0��ʾĿ��ֵ��ͬ��Ϊ��ͬ�⣬1��ʾ���߱�����ͬ������ͬ��
%
% ע�⣺
%   1.A0������������ĳ�ʼ��Ⱥ��
%   2.�����С�����⡣
%   3.֧�����У�permutation��������2-opt��inversion mutation����֧�ֶ����ƣ�binary��������1λ��ת��one bit flip��
%   4.��ͼֻ��������Ŀ�����Ŀ�꣬�Ұ�����ÿ����������������ʱ����
% 
% 2023.11.29  �����߽�Ļ��ƣ���Ӧ��ʹ��ģ��Ϊ1������������ӵõ��߽磬�����ӵõ��ı߽��������������ļнǲ����
% 2023.11.30  1. ��������ѡ��Al�����ŵĽ���Ϊ�����⵱ǰLS�Ķ���Ӧѡ�������⺯��ֵ��õģ���ԭ����CV_bestû�и���
%             2. ���ж����ʼ��ʱ���߽�������Ľ����ʼ������������������⺯��ֵ�ģ��ڽ��ڵĽ����ױ�������½�֧�䣬����update�����У���������н������Ľⲻ��������£�
%                Ҫ���������ұȵ�ǰ��������⺯��ֵ���õĽ��Ƿǳ��ѵģ���˻��ǽ���ʼ�⻮������
% 2023.12.1   ׼��1��2�����¼�����Ŀ�꺯��ֵ�������ظ����㣬������Nei�����������������

PROBLEM.Current(problem);
A0_new = [];
n = length(A0(1).dec);  % ������������
m = length(A0(1).obj);  % Ŀ�����
cheRefPoint = max(A0.objs,[],1);  % �ο���
[W, L] = UniformPoint(L, m);  % ����Ȩ�����������������������
log = [];
timecost = 0;

% ACN���
L_ACN = 10;
% drho = 0.1;

% ��������ӻ�
if process_visualization
    figure
    for i = 1:size(W,1)
        if i < size(W,1)
            subp_boundary = W(i,:)/norm(W(i,:)) + W(i+1,:)/norm(W(i+1,:));
            tmp = cheRefPoint - subp_boundary*min(cheRefPoint./subp_boundary);
            if m == 2
                hold on, plot([tmp(1) cheRefPoint(1)], [tmp(2) cheRefPoint(2)], 'k--');
            elseif m == 3
                view(3)
                grid on
                grid minor
                hold on, plot3([tmp(1) cheRefPoint(1)], [tmp(2) cheRefPoint(2)], [tmp(3) cheRefPoint(3)], 'k--');
            end
        end
        tmp = cheRefPoint - W(i,:)*min(cheRefPoint./W(i,:));
        if m == 2
            hold on, plot([tmp(1) cheRefPoint(1)], [tmp(2) cheRefPoint(2)],'k');
        elseif m == 3
            hold on, plot3([tmp(1) cheRefPoint(1)], [tmp(2) cheRefPoint(2)], [tmp(3) cheRefPoint(3)], 'k');
        end
    end
    title('original solution (red), every accepted solution (blue)')
    pause(0.01);
    clear tmp
end


% tmp1 = zeros(L*(L-1), m);  % ***����ҵ�����������
% tmp2 = zeros(L*(L-1), m);
% for i = 1 : L
%     tmp1((i-1)*(L-1)+1 : i*(L-1), :) = repmat(W(i,:),L-1,1);
%     tmp2((i-1)*(L-1)+1 : i*(L-1), :) = W([1:i-1 i+1:end], :);
% end
% cos = calCos(tmp1,tmp2);

for l = 1 : L
    start_time = tic;
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
            log = [log  ['��ʼʱ��������'  char(string(l))  '��û�н�'  newline]];
            clc
            disp(log);
            refPointMinusFit = repmat(cheRefPoint,length(A0),1) - A0.objs;
            [~, I] = sort(calCos(refPointMinusFit, repmat(myWV,length(A0),1)), 'descend');
            Al0 = A0(I == 1);
        end
    end
    
    % �������ʼ����ӻ�
    if process_visualization
        tmp = Al0.objs;
        if m == 2
            hold on, plot(tmp(:,1),tmp(:,2),'ro');
        elseif m == 3
            hold on, plot3(tmp(:,1),tmp(:,2),tmp(:,3),'ro');
        end
        pause(0.01);
        clear tmp
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
                % weiFitness = sum((cheRefPoint-Al(i_itSol).obj).*myWV, 2);  % ��Ϊ����󻯣����������ǳ˺�
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
                                update_plot_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                acceptCriIndex,solPrime,cheRefPoint,myWV,myWV_,rho,otherWV,strict_same,stepIndex,process_visualization, m);  % in
                            
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
                                    update_plot_ACN(inBoundCount,acceptFlag,Al0,explored,breakNeiExpFlag, ...  % in and out
                                    acceptCriIndex,solPrime,cheRefPoint,myWV,myWV_,rho,otherWV,strict_same,stepIndex,process_visualization, m);  % in
                                
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

                        % log_rho = [log_rho  rho_cur];
                        rho_cur
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
    region_timecost = toc(start_time);
    timecost = timecost + region_timecost;
    log = [log  ['������'  char(string(l))  '������'  str  '����ʱ'  char(string(region_timecost))  's'  newline]];
    clc
    disp(log);
    A0_new = [A0_new Al0];
end
log = [log  ['��������������ʱ'  char(string(timecost))  's'  newline]];
clc
disp(log);
A0_new = A0_new(NDSort(A0_new.objs, A0_new.cons, inf) == 1);  % ������з�֧������
end
