clear, close all


%% 参数设置
maxRuntime = 100;  % 子区域最大运行时间
L = 6;  % 子区域数量
n_randInitSol = 1;  % 随机初始化解的数目；1为原文设置
process_visualization = 1;  % 过程可视化，仅对CPLSD有用
strict_same = 0;  % 接受解时如何判断解相等，0目标向量相等，1决策向量相等（严格相等）；对于TSP，0为原文设置
is_recheck = 0;  % 0 1


%% 路径设置
pathstr = fileparts(mfilename('fullpath'));  % pathstr才是本m文件所在的路径
cd(pathstr);  % 更改当前活动目录路径
addpath(genpath(pathstr));
addpath(genpath('PlatEMO-master\PlatEMO'));


%% 问题设置
% binary
% pro = MONRP('N',n_randInitSol,'M',2,'D',100);  % 量纲不一致
% pro = MONRP2('N',n_randInitSol,'M',2,'D',100);  % 量纲一致
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


%% 初始解设置
% 随机初始化
rng(1)
PROBLEM.Current(pro);
A0 = pro.Initialization();

% 通过已保存数据初始化
% load 'PlatEMO格式的结果数据.mat'
% A0 = result{end};
% % A0 = UpdateArchive4paper1(A0,[],pro);
% W = UniformPoint(L, 3);
% [~,I] = min(pdist2(A0.objs-min(A0.objs,[],1),W,'cosine'),[],1);
% A0 = A0(unique(I));


%% 优化
% PPLS/D-C&A
K = 2;
% tp = 0.2:0.2:0.8;
tp = [1/3 2/3];
L_ACN = 100;
rho_init = 0;  % 为0时，初始为TCH；为1时，初始为WS
drho = 0;  % 为0时，关闭ACN
type_restart = 1;
type_c1 = 0;
[A0_new, norm_scope, log_rho, n_restart, log_tree] = CPPLSD_CA2_ACN(pro, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart, type_c1);

% PPLS/D-A    有bug未修复，详见CPPLSD_CA2_ACN中 2023.12.12 log
% K = 2;
% tp = [1/3 2/3];
% L_ACN = 30;
% rho_init = 1;
% drho = 0.1;
% type_restart = 4;  % 3感觉也不见得效果比2更好（实验条件：10s,L=6,NRP2 100 或 KP 250）
% [A0_new, norm_scope, log_rho, n_restart] = CPPLSD_CA_ACN(pro, A0, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart);

% PPLS/D
% A0_new = CPLSD_ACN(pro, A0, L, maxRuntime, process_visualization, strict_same, is_recheck, 0, 0);  % 可视化轨迹，不支持并行，注意画图消耗时间
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 1, 0.1);  % 并行版本  自适应
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 1, 0);  % 并行版本  WS
% A0_new = CPPLSD_ACN(pro, A0, L, maxRuntime, strict_same, is_recheck, 0, 0);  % 并行版本  TCH


%% 绘图对比
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
    % legend('LS后', 'LS前', '归一化', 'Location', 'best');

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
    legend('LS后', 'LS前', 'Location', 'best');
end

%% 计算HV
I_H = HV(A0_new,pro.optimum)

