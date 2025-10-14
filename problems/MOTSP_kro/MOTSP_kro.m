classdef MOTSP_kro < PROBLEM
% <multi/many> <permutation> <large/none>
% The multi-objective traveling salesman problem
% dataNo --- 1 --- 

%------------------------------- Reference --------------------------------
% TSPLIB & PPLS/D
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        C;  % Adjacency matrix of each map
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
            dataNo = obj.ParameterSet(1);  % 1: AB; 2: ABC; 3: ABCDE
            kroData;
            switch dataNo
                case 1
                    obj.M = 2;
                    dist_kroA = pdist2(kroA(:,2:3),kroA(:,2:3));
                    dist_kroB = pdist2(kroB(:,2:3),kroB(:,2:3));
                    obj.C{1} = dist_kroA;
                    obj.C{2} = dist_kroB;
                case 2
                    obj.M = 3;
                    dist_kroA = pdist2(kroA(:,2:3),kroA(:,2:3));
                    dist_kroB = pdist2(kroB(:,2:3),kroB(:,2:3));
                    dist_kroC = pdist2(kroC(:,2:3),kroC(:,2:3));
                    obj.C{1} = dist_kroA;
                    obj.C{2} = dist_kroB;
                    obj.C{3} = dist_kroC;
                case 3
                    obj.M = 5;
                    dist_kroA = pdist2(kroA(:,2:3),kroA(:,2:3));
                    dist_kroB = pdist2(kroB(:,2:3),kroB(:,2:3));
                    dist_kroC = pdist2(kroC(:,2:3),kroC(:,2:3));
                    dist_kroD = pdist2(kroD(:,2:3),kroD(:,2:3));
                    dist_kroE = pdist2(kroE(:,2:3),kroE(:,2:3));
                    obj.C{1} = dist_kroA;
                    obj.C{2} = dist_kroB;
                    obj.C{3} = dist_kroC;
                    obj.C{4} = dist_kroD;
                    obj.C{5} = dist_kroE;
                case 4
                    obj.M = 4;
                    dist_kroA = pdist2(kroA(:,2:3),kroA(:,2:3));
                    dist_kroB = pdist2(kroB(:,2:3),kroB(:,2:3));
                    dist_kroC = pdist2(kroC(:,2:3),kroC(:,2:3));
                    dist_kroD = pdist2(kroD(:,2:3),kroD(:,2:3));
                    obj.C{1} = dist_kroA;
                    obj.C{2} = dist_kroB;
                    obj.C{3} = dist_kroC;
                    obj.C{4} = dist_kroD;
            end
            obj.D = 100;
            obj.encoding = 'permutation';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [N,D]  = size(PopDec);
            PopObj = zeros(N,obj.M);
            for i = 1 : obj.M
                for j = 1 : N
                    for k = 1 : D-1
                        PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,k),PopDec(j,k+1));
                    end
                    PopObj(j,i) = PopObj(j,i) + obj.C{i}(PopDec(j,D),PopDec(j,1));
                end
            end
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = zeros(1,obj.M);
            for i = 1 : length(obj.C)
                R(i) = sum(max(obj.C{i},[],2));
            end
        end
    end
end