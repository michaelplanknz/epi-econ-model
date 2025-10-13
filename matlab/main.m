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

% Get model parameters
par = getPar();




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
[fu, results_u] = objFnCentFull(xu, par);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve centralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve central planner's problem with heuristic control function
x0 = [0, 0]; 
xOptHeur = fmincon(@(x)objFnCentHeur(x, par), x0, [], [], [], [], zeros(size(x0)), [1 inf] );

% Get model solution and costs for heuristic optimum
[f, resultsCentHeur] = objFnCentHeur(xOptHeur, par);

% Use heuristic solution to set initial condition for full optimization
x0 = interp1(t', resultsCentHeur.a', tMesh)';

% Solve central planner's problem with full time-dependent control variable
xOptCent = fmincon(@(x)objFnCentFull(x, par), x0, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups) );

% Get model solution and costs for full optimum
[~, resultsCent] = objFnCentFull(xOptCent, par);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve decentralised problem with analytic method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise at uncontrolled state
a0 = ones(par.nGroups, length(t));
resultsDecent0 = getResultsAnalytic("decent", a0, par);


% Initialise at centralised solution
a0 = resultsCent.a;
resultsDecent1 = getResultsAnalytic("decent", a0, par);


objFnDecentFull(resultsDecent0.a, resultsDecent1.a, par)
objFnDecentFull(resultsDecent1.a, resultsDecent0.a, par)


fOut = outFolder + "results.mat";
save(fOut);


