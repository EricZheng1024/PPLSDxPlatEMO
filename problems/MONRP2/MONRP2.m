classdef MONRP2 < PROBLEM
% <multi> <binary> <large/none>
% The multi-objective next release problem with normalization
% m --- 100 --- Number of customers

%------------------------------- Reference --------------------------------
% Y. Zhang, M. Harman, and S. A. Mansouri, The multi-objective next release
% problem, Proceedings of the 9th Annual Conference on Genetic and
% Evolutionary Computation, 2007, 1129-1137.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        Cost;   % Cost of each requirement
        Value;  % Value of each customer on each requirement
        sumCost;
        sumValue;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            m = obj.ParameterSet(100);
            obj.M = 2;
            if isempty(obj.D); obj.D = 100; end
            obj.encoding = 'binary';
            % Randomly generate costs and values
            n    = obj.D;
            file = sprintf('MONRP2-n%d-m%d.mat',n,m);
            file = fullfile(fileparts(mfilename('fullpath')),file);
            if exist(file,'file') == 2
                load(file,'Cost','Value');
            else
                Cost  = randi(9,1,n);
                Value = randi([0 5],n,m);
                save(file,'Cost','Value');
            end
            obj.Cost  = Cost;
            obj.Value = Value;
            obj.sumCost = sum(obj.Cost);
            obj.sumValue = sum(obj.Value(:));
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj(:,1) = sum(repmat(obj.Cost,size(PopDec,1),1).*PopDec,2) / obj.sumCost;
            PopObj(:,2) = (obj.sumValue - sum(PopDec*obj.Value,2)) / obj.sumValue;
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(~,~)
            R = [1,1];
        end
        %% Display a population in the objective space
        function DrawObj(obj,Population)
            PopObj(:,1) = sum(repmat(obj.Cost,length(Population),1).*Population.decs,2) / obj.sumCost;
            PopObj(:,2) = (obj.sumValue - sum(Population.decs*obj.Value,2)) / obj.sumValue;
            Draw(PopObj,{'Cost(normalization)','Dissatisfaction score(normalization)',[]});
        end
    end
end