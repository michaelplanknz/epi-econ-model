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

% % Solve central planner's problem with heuristic control function
% x0 = xu; 
% xOptCent = fmincon(@(x)objFnCent(x, par), x0, [], [], [], [], zeros(size(x0)), [] );
% 
% % Get model solution and costs for heuristic optimum
% [f, resultsCent] = objFnCent(xOptCent, par);

% Solve analytic approximation 
resultsCentAnalytic = getResultsAnalytic("cent", par);

% Set initial condition for the full optimization by interpolating the
% analytic approximatoin at the time mesh points

x0 = interp1(t', resultsCentAnalytic.a', tMesh)';

% Solve central planner's problem with full time-dependent control variable
xOptCent = fmincon(@(x)objFnCentFull(x, par), x0, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups) );

% Get model solution and costs for full optimum
[~, resultsCent] = objFnCentFull(xOptCent, par);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve decentralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Start by setting background behaviour to the solution to the centralised
% % problem with heuristic control function
% xOptDecent = xOptCent;
% 
% iRep = 1;
% convFlag = false;
% while iRep <= maxReps & ~convFlag 
%     xSav = xOptDecent;
%     % Optimize each group in turn 
%     for iGroup = 1:par.nGroups
%         % Indices of x relating to the current group
%         groupInd = [1, 2] + 2*(iGroup-1);
%         xOptGroup = fmincon(@(x)objFnDecent(x, xOptDecent, iGroup, par), xOptDecent(groupInd), [], [], [], [], zeros(2, 1), [] );
%         xOptDecent(groupInd) = xOptGroup;
%     end
%     convFlag = norm(xOptDecent-xSav)/norm(xSav) < relTol;
%     iRep = iRep+1;
% end
% 
% if ~convFlag
%     fprintf('Warning: decentralised problem with heuristic control function has failed to converge after %i iterations\n', maxReps)
% end
% 
% % Get model solution and costs for heuristic optimum
% [f, resultsDecent] = objFnDecentBoth(xOptDecent, xOptDecent, par);
% 
% % Set initial condition for the full optimization by interpolating the
% % heuristic function and the time mesh points
% xOptDecentFull = interp1(t', resultsDecent.a', tMesh)';
% 
% iRep = 1;
% convFlag = false;
% while iRep <= maxReps & ~convFlag 
%     % Start by setting background behaviour to the solution central problem
%     xSav = xOptDecentFull;
%     % Optimize each group in turn 
%     for iGroup = 1:par.nGroups
%         % Indices of x relating to the current group
%         xOptGroup = fmincon(@(x)objFnDecentFull(x, xOptDecentFull, iGroup, par), xOptDecentFull(:, iGroup), [], [], [], [], zeros(nMeshPts, 1), ones(nMeshPts, 1) );
%         xOptDecentFull(:, iGroup) = xOptGroup;
%     end
%     convFlag = norm(xOptDecentFull-xSav)/norm(xSav) < relTol;
%     iRep = iRep+1;
% end
% 
% if ~convFlag
%     fprintf('Warning: decentralised problem with full control function has failed to converge after %i iterations\n', maxReps)
% end
% 
% xOptDecentFull = xOptDecentFull';
% 
% % Get model solution and costs for full optimum
% [~, resultsDecentFull] = objFnDecentFullBoth(xOptDecentFull, xOptDecentFull, par);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve centralised and decentralised problem with Shaun's method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsDecent = getResultsAnalytic("decent", par);


fOut = outFolder + "results.mat";
save(fOut);


