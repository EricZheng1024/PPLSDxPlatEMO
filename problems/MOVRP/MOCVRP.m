classdef MOCVRP < PROBLEM
% <multi> <permutation> <large/none> <constrained>
% The multi-objective vehicle routing problem
% 2023.6.15 订正了保存的数据文件名
% c --- 0 --- Correlation parameter
% Nc --- 30 --- the number of clients(excluding warehouse)
% Nv --- 10 --- the number of vehicle
% Dm --- 0.2 --- the max demand of clients(capacity of vehicle is 1)

% 带容量约束的多车VRP
% 设仓库在C和demand中为客户编号1，车辆编号在客户编号之后
% 决策变量索引1代表的为C和demand中的客户编号2

    properties(Access = private)
        C;  % Adjacency matrix of each map
        Nv;
        demand;  % the demand of clients
        v_no;  %  vehicle numbers
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            [c, Nc, obj.Nv, Dm] = obj.ParameterSet(0, 30, 10, 0.2);
            if isempty(obj.M); obj.M = 2; end
            obj.D = Nc + obj.Nv - 1;
            obj.encoding = 'permutation';
            % Randomly generate the adjacency matrices
            file = sprintf('MOCVRP-M%d-D%d-c%g-Nc%d-Nv%d-Dm%.2f.mat',obj.M,obj.D, c,Nc,obj.Nv,Dm);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'C','demand');
            else
                C = cell(1,obj.M);
                C{1} = rand(Nc+1);  % 包含仓库
                for i = 2 : obj.M
                    C{i} = c*C{i-1} + (1-c)*rand(Nc+1);
                end
                for i = 1 : obj.M
                    C{i} = tril(C{i},-1) + triu(C{i}',1);
                end
                demand = [0, rand(1,Nc) * Dm];
                save(file,'C','demand');
            end
            obj.C = C;
            obj.demand = demand;
            obj.v_no = Nc+1 : Nc+obj.Nv;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            PopObj = zeros(N,obj.M);
            for i = 1 : obj.M
                for j = 1 : N
                    x = [0 PopDec(j,:)];
                    x(ismember(x, obj.v_no)) = 0;
                    x = x + 1;
                    PopObj(j,i) = trace(obj.C{i}(x(1:end-1), x(2:end)));
                end
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
                    PopCon(i,j) = sum(obj.demand(X(i,index(j)+1:index(j+1)-1)+1), 2) - 1;  % 车辆容量为1
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