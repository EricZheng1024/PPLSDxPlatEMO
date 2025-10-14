clear, close all


%% ��������
maxRuntime = 100;  % �������������ʱ��
L = 6;  % ����������
n_randInitSol = 1;  % �����ʼ�������Ŀ��1Ϊԭ������
process_visualization = 1;  % ���̿��ӻ�������CPLSD����
strict_same = 0;  % ���ܽ�ʱ����жϽ���ȣ�0Ŀ��������ȣ�1����������ȣ��ϸ���ȣ�������TSP��0Ϊԭ������
is_recheck = 0;  % 0 1


%% ·������
pathstr = fileparts(mfilename('fullpath'));  % pathstr���Ǳ�m�ļ����ڵ�·��
cd(pathstr);  % ���ĵ�ǰ�Ŀ¼·��
addpath(genpath(pathstr));
addpath(genpath('PlatEMO-master\PlatEMO'));


%% ��������
% binary
% pro = MONRP('N',n_randInitSol,'M',2,'D',100);  % ���ٲ�һ��
% pro = MONRP2('N',n_randInitSol,'M',2,'D',100);  % ����һ��
% pro = MOKP('N',n_randInitSol,'M',2,'D',250);
% pro = MOKP('N',n_randInitSol,'M',2,'D',1000);

% permutation
% pro = MOTSP('N',n_randInitSol,'M',2,'D',30,'parameter',{0});
% pro = MOTSP('N',n_randInitSol,'M',3,'D',100);
pro = MOTSP_kro('N',n_randInitSol,'M',2,'D',100,'parameter',{1});
% pro = MOTSP_other('N',n_randInitSol,'parameter',{1});

% pro = mQAP('N',n_randInitSol,'M',2,'D',10);
% pro = mQAP('N',n_randInitSol,'M',3,'D',100);

% pro = MOCVRP2('N',n_randInitSol,'M',2,'parameter',{20,4,5,20});
% pro = MOCVRP('N',n_randInitSol,'M',2,'parameter',{0,20,10,0.2});
% pro = MOCVRP('N',n_randInitSol,'M',3,'parameter',{0,20,10,0.2});
% pro = MOCVRP('N',n_randInitSol,'M',2,'parameter',{0,10,5,0.4});
% pro = MOCVRP('N',n_randInitSol,'M',3,'parameter',{0,10,5,0.4});


%% ��ʼ������
% �����ʼ��
rng(1)
PROBLEM.Current(pro);
A0 = pro.Initialization();

% ͨ���ѱ������ݳ�ʼ��
% load 'PlatEMO��ʽ�Ľ������.mat'
% A0 = result{end};
% % A0 = UpdateArchive4paper1(A0,[],pro);
% W = UniformPoint(L, 3);
% [~,I] = min(pdist2(A0.objs-min(A0.objs,[],1),W,'cosine'),[],1);
% A0 = A0(unique(I));


%% �Ż�
% PPLS/D-C&A
K = 2;
% tp = 0.2:0.2:0.8;
tp = [1/3 2/3];
L_ACN = 100;
rho_init = 0;  % Ϊ0ʱ����ʼΪTCH��Ϊ1ʱ����ʼΪWS
drho = 0;  % Ϊ0ʱ���ر�ACN
type_restart = 1;
type_c1 = 0;
[A0_new, norm_scope, log_rho, n_restart, log_tree] = CPPLSD_CA2_ACN(pro, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart, type_c1);

% PPLS/D-A    ��bugδ�޸������CPPLSD_CA2_ACN�� 2023.12.12 log
% K = 2;
% tp = [1/3 2/3];
% L_ACN = 30;
% rho_init = 1;
% drho = 0.1;
% type_restart = 4;  % 3�о�Ҳ������Ч����2���ã�ʵ��������10s,L=6,NRP2 100 �� KP 250��
% [A0_new, norm_scope, log_rho, n_restart] = CPPLSD_CA_ACN(pro, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart);

% PPLS/D
% A0_new = CPLSD_ACN(pro, A0, L, maxRuntime, process_visualization, strict_same, is_recheck, 0, 0);  % ���ӻ��켣����֧�ֲ��У�ע�⻭ͼ����ʱ��
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 1, 0.1);  % ���а汾  ����Ӧ
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 1, 0);  % ���а汾  WS
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 0, 0);  % ���а汾  TCH


%% ��ͼ�Ա�
close all
figure
after = A0_new.objs;
before = A0.objs;
switch pro.M
    case 2
        plot(after(:, 1), after(:, 2), 'b*');
        hold on, plot(before(:, 1), before(:, 2), 'ro');
        
    case 3
        plot3(after(:, 1), after(:, 2), after(:, 3), 'b*');
        hold on, plot3(before(:, 1), before(:, 2), before(:, 3), 'ro');
        grid on
end

try
    % plot(norm_scope(:,1),norm_scope(:,2), 's')
    % legend('LS��', 'LSǰ', '��һ��', 'Location', 'best');

    % figure
    % hold on
    % plot(log_rho{1})
    % plot(log_rho{round(length(log_rho)/2)})
    % plot(log_rho{end})
    % legend

    figure
    hold on
    for i = 1 : length(log_tree)
        objs = log_tree{i}.objs;
        plot(objs(:,1), objs(:,2), '*')
    end
catch
    legend('LS��', 'LSǰ', 'Location', 'best');
end

%% ����HV
I_H = HV(A0_new,pro.optimum)

