clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set output folder location
outFolder = '../output/';

% Numerical parameters
% Maximum number of iterates and convergence tolerance for calculating Nash equilibrium
maxReps = 100;
relTol = 1e-4;

% Get default model parameters
par = getPar();

% Define parameters to vary
Beta_vals = 0.25:0.05:0.7;
costPerInf_vals = 0.1:0.1:1.5;

% Get a list of parameter combinations
[Beta_mat, costPerInf_mat] = meshgrid(Beta_vals, costPerInf_vals);
Beta_list = Beta_mat(:);
costPerInf_list = costPerInf_mat(:);
nScenarios = length(Beta_list);

% Set options for opimization routine 
%opts1 = optimoptions('fmincon', 'Display', 'notify');
opts = optimoptions('fmincon', 'Display', 'notify', 'MaxFunctionEvaluations', 20000, 'MaxIterations', 5000);

% Loop through parameter combinations
for iScenario = 1:nScenarios

    fprintf('Scenario %i/%i\n', iScenario, nScenarios)

    %Set scneario-dependent parameters
    par.Beta = Beta_list(iScenario);
    par.costPerInf = costPerInf_list(iScenario);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve unmitigated epidemic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set up time vector, and time mesh points for defining control function
    t = 0:par.dt:par.tMax;
    tMesh = 0:par.meshSpace:par.tMax;

    nTimePts = length(t);
    
    % Define no control state
    nMeshPts = length(tMesh);
    xu = ones(par.nGroups, nMeshPts); 
    
    % Get model solution and costs for no control 
    [~, results_u(iScenario)] = objFnCentFull(xu, par);
    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial conditions for centralised and decentralised
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if iScenario == 1
        % Need to initialise optimization routine somehow on first
        % iteration

        % Initialise centralised problem at uncontrolled state
        x0 = ones(1, length(tMesh));

        % Initialise decentralised problem at uncontrolled state
        a0 = ones(par.nGroups, length(t));
        a0_SD = ones(par.nGroups, length(t));
    else
        % On subsequent iterations use a previous solution to initialise
        if Beta_list(iScenario) == Beta_list(iScenario-1)
            % Use the most recdent solution if it had the same value of Beta
            ind = iScenario-1;
        else
            % Otherwise use the most recent solution that had the same
            % value of costPerInf
            ind = find(costPerInf_list(1:iScenario-1) == costPerInf_list(iScenario), 1, 'last');
        end
        x0 = resultsCent(ind).x;
        a0 = resultsDecent(ind).a;
        a0_SD = resultsDecent_SD(ind).a;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve centralised problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve central planner's problem with full time-dependent control variable
    xOptCent = fmincon(@(x)objFnCentFull(x, par), x0, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups), [], opts );
    
    % Get model solution and costs for full optimum
    [~, resultsCent(iScenario)] = objFnCentFull(xOptCent, par);
    
    % Use longer time period if solution end up under HIT
    R0 = eigs(par.Beta, [], 1)/par.Gamma;
    if resultsCent(iScenario).S(end)/par.N > 1/R0
        fprintf('    Warning: final size %.3f is under HIT %.3f, may need to run for longer\n', 1-resultsCent(iScenario).S(end)/par.N, 1-1/R0)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve decentralised problem with analytic method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    resultsDecent(iScenario) = getResultsAnalytic(a0, par);

    resultsDecent_SD(iScenario) = getResultsAnalytic_SD(a0_SD, par);

end

fOut = outFolder + "results.mat";
save(fOut);


