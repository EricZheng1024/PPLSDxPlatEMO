classdef MOCVRP2 < PROBLEM
% <multi> <permutation> <large/none> <constrained>
% MOCVRP，欧拉图C，目标1为路径总长最小，目标2为车辆的最长路径最短
% 默认参数：50+1个节点，每个节点的demand为1~9的随机数（depot除外）
% 设depot在C和demand中为客户编号1，车辆编号在客户编号之后
% 决策变量索引1代表的为C和demand中的客户编号2
% 2023.6.15 订正了保存的数据文件名
% Nc --- 50 --- the number of clients (excluding depot)
% Dm --- 9 --- the max demand of clients (1~Dm integer)
% Nv --- 12 --- the number of vehicle
% capa --- 40 --- the capacity of vehicle

    properties(Access = private)
        C;  % Adjacency matrix of each map
        demand;  % Demand of clients
        Nv;  % Vehicle 数量
        v_no;  %  Vehicle 编号
        capa;  % Vehicle capacity
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [Nc, Dm, obj.Nv, obj.capa] = obj.ParameterSet(50, 9, 12, 40);
            obj.M = 2;
            obj.D = Nc + obj.Nv - 1;
            obj.encoding = 'permutation';
            % Randomly generate the adjacency matrices
            file = sprintf('MOCVRP2-M%d-D%d-Nc%d-Dm%d-Nv%d-capa%d.mat',obj.M,obj.D, Nc,Dm,obj.Nv,obj.capa);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'C','demand');
            else
                map = rand(Nc+1,2);  % 包含仓库
                C = pdist2(map, map);
                demand = [0 randi(Dm,1,Nc)];  % 1~Dm, 0 for depot
                save(file,'map','C','demand');
            end
            obj.C = C;
            obj.demand = demand;
            obj.v_no = Nc+1 : Nc+obj.Nv;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,~]  = size(PopDec);
            PopObj = zeros(N,obj.M);
            for i = 1 : N
                x = [0 PopDec(i,:) 0];
                x(ismember(x, obj.v_no)) = 0;
                x = x + 1;
%                 PopObj(i,1) = trace(obj.C(x(1:end-1), x(2:end)));  % 总距离
                v_index = find(x == 1);
                tmp = zeros(1,length(v_index)-1);
                for j = 1 : length(v_index)-1
                    tmp(j) = trace(obj.C(x(v_index(j):(v_index(j+1)-1)), x((v_index(j)+1):v_index(j+1))));
                end
                PopObj(i,1) = sum(tmp);  % 总距离（与上面等价）
                PopObj(i,2) = max(tmp);  % 车辆最长路径距离
            end
        end
        %% Calculate constraint violations
        function PopCon = CalCon(obj,X)
            [N,D]  = size(X);
            X = [(D+1)*ones(N,1) X (D+1)*ones(N,1)];
            PopCon = zeros(N,obj.Nv);
            for i = 1 : N
                index = find(ismember(X(i,:), obj.v_no));
                for j = 1 : obj.Nv
                    PopCon(i,j) = sum(obj.demand(X(i,index(j)+1:index(j+1)-1)+1), 2) - obj.capa;
                end
            end
            PopCon(PopCon < 0) = 0;
            PopCon = sum(PopCon, 2);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = zeros(1,obj.M) + obj.D;
        end
    end
end