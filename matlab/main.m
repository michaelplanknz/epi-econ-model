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
Beta_arr       = [0.3   0.3   0.3   0.6   0.6   0.6];
costPerInf_arr = [0.1   0.5   1.2   0.1   0.5   1.2];
nScenarios = length(Beta_arr);

for iScenario = 1:nScenarios

    %Set scneario-dependent parameters
    par.Beta = Beta_arr(iScenario);
    par.costPerInf = costPerInf_arr(iScenario);

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
    
    % OLD METHOD (needed if initialising at analytical approximation
    % doesn't work)
    % Solve central planner's problem with heuristic control function
    x0 = [0, 0]; 
    xOptHeur = fmincon(@(x)objFnCentHeur(x, par), x0, [], [], [], [], zeros(size(x0)), [1 inf] );
    
    % Get model solution and costs for heuristic optimum
    [f, resultsCentHeur(iScenario)] = objFnCentHeur(xOptHeur, par);
    
    % Use heuristic solution to set initial condition for full optimization
    x0 = interp1(t', resultsCentHeur(iScenario).a', tMesh)';
    
    % Initialise at uncontrolled state
%      a0 = ones(par.nGroups, length(t));
%      % Get analytical approximation to optimal centralised control
%     results_tmp = getResultsAnalytic("cent", a0, par);
%      % Use this to set initial condition for full optimization routine
%      x0 = interp1(t, results_tmp.a, tMesh);

    % Solve central planner's problem with full time-dependent control variable
    xOptCent = fmincon(@(x)objFnCentFull(x, par), x0, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups) );
    
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


