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
% Compute elimination costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate linear and quadratic costs per unit time if same control is
% applied to whole pop
costlin_avg = par.controlFrac * sum(par.costlin.*par.N);
costquad_avg = par.controlFrac * sum(par.costquad.*par.N);


% Calculate optimal a for elimination of small outbreaks
R0 = eigs(par.Beta, [], 1)/par.Gamma;
myF = @(a)( (costlin_avg*(1-a) + costquad_avg*(1-a).^2)./(par.Gamma - par.Gamma*R0*par.alpha_TTI*a.^2) );
aOptElim = fmincon(myF, 0, [], [], [], [], 0, sqrt(1/(par.alpha_TTI*R0)) );

% Calculate elimination cost per unit time
Celim_per_unit_time = par.b + 1/3 * par.r*log(par.xOutbreak)*(costlin_avg*(1-aOptElim) + costquad_avg*(1-aOptElim)^2)/(par.Gamma - R0*par.Gamma*par.alpha_TTI*aOptElim^2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve unmitigated epidemic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up time vector, and time mesh points for defining control function
t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

% Define no control state
xu = repmat([0; 0], par.nGroups, 1);

% Get model solution and costs for no control 
[fu, results_u] = objFnCent(xu, par);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve centralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve central planner's problem with heuristic control function
x0 = xu; 
xOptCent = fmincon(@(x)objFnCent(x, par), x0, [], [], [], [], zeros(size(x0)), [] );

% Get model solution and costs for heuristic optimum
[f, resultsCent] = objFnCent(xOptCent, par);

% Set initial condition for the full optimization by interpolating the
% heuristic function and the time mesh points
nMeshPts = length(tMesh);
x0Full = interp1(t', resultsCent.a', tMesh)';

% Solve central planner's problem with full time-dependent control variable
xOptCentFull = fmincon(@(x)objFnCentFull(x, par), x0Full, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups) );

% Get model solution and costs for full optimum
[~, resultsCentFull] = objFnCentFull(xOptCentFull, par);

% Group labels for legend
groupLbls = string((1:par.nGroups)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve decentralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start by setting background behaviour to the solution to the centralised
% problem with heuristic control function
xOptDecent = xOptCent;

iRep = 1;
convFlag = false;
while iRep <= maxReps & ~convFlag 
    xSav = xOptDecent;
    % Optimize each group in turn 
    for iGroup = 1:par.nGroups
        % Indices of x relating to the current group
        groupInd = [1, 2] + 2*(iGroup-1);
        xOptGroup = fmincon(@(x)objFnDecent(x, xOptDecent, iGroup, par), xOptDecent(groupInd), [], [], [], [], zeros(2, 1), [] );
        xOptDecent(groupInd) = xOptGroup;
    end
    convFlag = norm(xOptDecent-xSav)/norm(xSav) < relTol;
    iRep = iRep+1;
end

if ~convFlag
    fprintf('Warning: decentralised problem with heuristic control function has failed to converge after %i iterations\n', maxReps)
end

% Get model solution and costs for heuristic optimum
[f, resultsDecent] = objFnDecentBoth(xOptDecent, xOptDecent, par);

% Set initial condition for the full optimization by interpolating the
% heuristic function and the time mesh points
xOptDecentFull = interp1(t', resultsDecent.a', tMesh)';

iRep = 1;
convFlag = false;
while iRep <= maxReps & ~convFlag 
    % Start by setting background behaviour to the solution central problem
    xSav = xOptDecentFull;
    % Optimize each group in turn 
    for iGroup = 1:par.nGroups
        % Indices of x relating to the current group
        xOptGroup = fmincon(@(x)objFnDecentFull(x, xOptDecentFull, iGroup, par), xOptDecentFull(:, iGroup), [], [], [], [], zeros(nMeshPts, 1), ones(nMeshPts, 1) );
        xOptDecentFull(:, iGroup) = xOptGroup;
    end
    convFlag = norm(xOptDecentFull-xSav)/norm(xSav) < relTol;
%    norm(xOptDecentFull-xSav)/norm(xSav)
%    figure(10);
%    plot(tMesh, xOptDecentFull)
%    drawnow
    iRep = iRep+1;
end

if ~convFlag
    fprintf('Warning: decentralised problem with full control function has failed to converge after %i iterations\n', maxReps)
end

xOptDecentFull = xOptDecentFull';

% Get model solution and costs for full optimum
[~, resultsDecentFull] = objFnDecentFullBoth(xOptDecentFull, xOptDecentFull, par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve centralised and decentralised problem with Shaun's method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsCentShaun = getResultsShaun("cent", par);
resultsDecentShaun = getResultsShaun("decent", par);


fOut = outFolder + "results.mat";
save(fOut);


