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

% Set options for opimization routine (for heuristic and full problems)
opts1 = optimoptions('fmincon', 'Display', 'notify');
opts2 = optimoptions('fmincon', 'Display', 'notify', 'MaxFunctionEvaluations', 20000, 'MaxIterations', 5000);

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
    
    % Define no control state
    nMeshPts = length(tMesh);
    xu = ones(par.nGroups, nMeshPts); 
    
    % Get model solution and costs for no control 
    [~, results_u(iScenario)] = objFnCentFull(xu, par);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve centralised problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iScenario == 1
        % Need to initialise optimization routine somehow on first
        % iteration

        % Solve central planner's problem with heuristic control function
        x0 = [0, 0]; 
        xOptHeur = fmincon(@(x)objFnCentHeur(x, par), x0, [], [], [], [], zeros(size(x0)), [1 inf], [], opts1 );
        
        % Get model solution and costs for heuristic optimum
        [f, resultsCentHeur(iScenario)] = objFnCentHeur(xOptHeur, par);
        
        % Use heuristic solution to set initial condition for full optimization
        x0 = interp1(t', resultsCentHeur(iScenario).a', tMesh)';
    else
        % On subsequent iterations use a previous solution to initialise
        if Beta_list(iScenario) == Beta_list(iScenario-1)
            % Use the most recdent solution if it had the same value of Beta
            x0 = resultsCent(iScenario-1).x;
        else
            % Otherwise use the most recent solution that had the same
            % value of costPerInf
            ind = find(costPerInf_list(1:iScenario-1) == costPerInf_list(iScenario), 1, 'last');
            x0 = resultsCent(ind).x;
        end
    end

    % Solve central planner's problem with full time-dependent control variable
    xOptCent = fmincon(@(x)objFnCentFull(x, par), x0, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups), [], opts2 );
    
    % Get model solution and costs for full optimum
    [~, resultsCent(iScenario)] = objFnCentFull(xOptCent, par);
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve decentralised problem with analytic method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialise at uncontrolled state
    a0 = ones(par.nGroups, length(t));
    resultsDecent(iScenario) = getResultsAnalytic("decent", a0, par);

end

fOut = outFolder + "results.mat";
save(fOut);


