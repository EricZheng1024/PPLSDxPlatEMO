classdef CPPLSD_CA2_ACN_PlatEMOALG < ALGORITHM
% <multi/many> <binary/permutation> <constrained/none>
% PPLS/D-C&A
% N --- 1 --- 
% L ---  --- 
% K ---  --- 
% rho_init --- 1 --- 
% drho --- 0.1 --- 
% L_ACN --- 100 --- 
% tp --- [1/3 2/3] --- 
% maxRuntime --- 300 --- 
% type_restart --- 1 --- 

%------------------------------- Reference --------------------------------
% 
%------------------------------- Copyright --------------------------------
% 
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [N, L, K, rho_init, drho, L_ACN, tp, maxRuntime, type_restart] = Algorithm.ParameterSet(1, [], [], 1, 0.1, 100, [1/3 2/3], 300, 1);
            switch Problem.M
                case 2
                    L = 20;
                    K = 2;
                case 3
                    L = 45;
                    K = 6;
            end
            strict_same = 0;
            N_tmp = Problem.N;
            Problem.N = floor(N);
            Population = Problem.Initialization();
            Problem.N = N_tmp;
            if maxRuntime < 0
                switch maxRuntime
                case -1
                    switch class(Problem)
                        case 'MONRP2'
                            maxRuntime = 10;
                        case 'MOKP'
                            switch Problem.D
                                case 600
                                    maxRuntime = 5;
                                case 1000
                                    maxRuntime = 10;
                            end
                        case {'MOTSP', 'mQAP'}
                            switch Problem.D
                                case 100
                                    maxRuntime = 10;
                                case 200
                                    maxRuntime = 20;
                            end
                    end
                end
            end
            type_c1 = abs((N - floor(N)) - 0.1) < 1e-12;  % 0(false)为从准则1*，0.1(true)为准则1

            %% Optimization
            while Algorithm.NotTerminated(Population)
                [Population, norm_scope, log_rho, n_restart] = CPPLSD_CA2_ACN_2022a(Problem, Population, L, K, rho_init, drho, L_ACN, tp, maxRuntime, strict_same, type_restart, type_c1);
                Population(end).appendix = {norm_scope, log_rho, n_restart};  % 使用最后一个解来保存信息
            end
        end
    end
end