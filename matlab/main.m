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
resultsDecent1 = getResultsAnalytic("decent", a0, par);


% Initialise at centralised solution
a0 = resultsCent.a;
resultsDecent2 = getResultsAnalytic("decent", a0, par);

a1 = resultsDecent1.a;
a2 = resultsDecent2.a;

w11 = objFnDecentFull(a1, a1, par)
w21 = objFnDecentFull(a2, a1, par)
w22 = objFnDecentFull(a2, a2, par)
w12 = objFnDecentFull(a1, a2, par)


z = 0:0.01:1.5;

for ii = 1:length(z)
    wz1(ii) = objFnDecentFull(1 - z(ii)*(1-a1), a1, par);
    wz2(ii) = objFnDecentFull(1 - z(ii)*(1-a2), a2, par);
end

figure;
plot(z, wz1, z, wz2);



fOut = outFolder + "results.mat";
save(fOut);


