clear 
close all

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
[fu, costInfu, costContu, ~, Su, Iu] = objFnCent(xu, par);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve centralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve central planner's problem with heuristic control function
x0 = xu; 
xOpt = fmincon(@(x)objFnCent(x, par), x0, [], [], [], [], zeros(size(x0)), [] );

% Get model solution and costs for heuristic optimum
[f, costInf, costCont, aOpt, S, I] = objFnCent(xOpt, par);

% Set initial condition for the full optimization by interpolating the
% heuristic function and the time mesh points
nMeshPts = length(tMesh);
x0Full = interp1(t', aOpt', tMesh)';

% Solve central planner's problem with full time-dependent control variable
xOptFull = fmincon(@(x)objFnCentFull(x, par), x0Full, [], [], [], [], zeros(nMeshPts, par.nGroups), ones(nMeshPts, par.nGroups) );

% Get model solution and costs for full optimum
[~, costInfFull, costContFull, aFull, SFull, IFull] = objFnCentFull(xOptFull, par);

% Calculate recovered compartment
R = par.N-S-I;
Ru = par.N-Su-Iu;
RFull = par.N-SFull-IFull;

groupLbls = string((1:par.nGroups)');

h = figure(1);
h.Position = [     240         125        1024         744];
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
plot(t, RFull)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, Ru, '--')
plot(t, aFull, 'LineWidth', 2)
if par.nGroups > 1
    lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
else
    lbls = ["R", "R_u", "a"];
end
legend(lbls, 'Location', 'southeast')
grid on
xlabel('time (days)')

nexttile;
plot(t, costInfFull+costContFull)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, costInfu+costContu, '--')
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["cost", "cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

nexttile;
plot(t, costInfFull, t, costContFull)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, costInfu, '--')
if par.nGroups > 1
    lbls = ["inf cost_"+groupIDs; "cont cost_"+groupIDs; "inf cost_"+groupIDs+"_u"];
else
    lbls = ["inf cost", "cont cost", "inf cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

sgtitle('centralized control')
drawnow




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve decentralised problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iRep = 1;
convFlag = false;
while iRep <= maxReps & ~convFlag 
    % Start by setting background behaviour to the solution central problem
    xSav = xOpt;
    % Optimize each group in turn 
    for iGroup = 1:par.nGroups
        % Indices of x relating to the current group
        groupInd = [1, 2] + 2*(iGroup-1);
        xOptGroup = fmincon(@(x)objFnDecent(x, xOpt, iGroup, par), xOpt(groupInd), [], [], [], [], zeros(2, 1), [] );
        xOpt(groupInd) = xOptGroup;
    end
    convFlag = norm(xOpt-xSav)/norm(xSav) < relTol;
    iRep = iRep+1;
end

% Get model solution and costs for heuristic optimum
[f, costInf_pp, costCont_pp, aOpt, S, I, S_focal] = objFnDecentBoth(xOpt, xOpt, par);

% Set initial condition for the full optimization by interpolating the
% heuristic function and the time mesh points
xOptFull = interp1(t', aOpt', tMesh)';

iRep = 1;
convFlag = false;
while iRep <= maxReps & ~convFlag 
    % Start by setting background behaviour to the solution central problem
    xSav = xOptFull;
    % Optimize each group in turn 
    for iGroup = 1:par.nGroups
        % Indices of x relating to the current group
        xOptGroup = fmincon(@(x)objFnDecentFull(x, xOptFull, iGroup, par), xOptFull(:, iGroup), [], [], [], [], zeros(nMeshPts, 1), ones(nMeshPts, 1) );
        xOptFull(:, iGroup) = xOptGroup;
    end
    convFlag = norm(xOptFull-xSav)/norm(xSav) < relTol;
%    norm(xOptFull-xSav)/norm(xSav)
%    figure(10);
%    plot(tMesh, xOptFull)
%    drawnow
    iRep = iRep+1;
end

xOptFull = xOptFull';

% Get model solution and costs for full optimum
[costInfFullDecent, costContFullDecent, aFullDecent, SFullDecent, IFullDecent] = deal(zeros(par.nGroups, length(t)));
for iGroup = 1:par.nGroups
    [~, costInfFullDecent(iGroup, :), costContFullDecent(iGroup, :), aFullDecent(iGroup, :), SFullDecent, IFullDecent] = objFnDecentFull(xOptFull(iGroup, :), xOptFull, iGroup, par);
end

% Calculate recovered compartment
R = par.N-S-I;
RFullDecent = par.N-SFullDecent-IFullDecent;

% Plotting
h = figure(2);
h.Position = [    360         183        1024         755];
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
plot(t, RFullDecent)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, Ru, '--')
plot(t, aFullDecent, 'LineWidth', 2)
grid on
legend('R', 'R_u', 'a', 'Location', 'southeast')
xlabel('time (days)')


nexttile;
plot(t, (costInfFullDecent+costContFullDecent).*par.N)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, costInfu+costContu, '--')
grid on
legend('cost', 'cost_u', 'Location', 'southeast')
ylabel('cumulative aggregate cost')
xlabel('time (days)')


nexttile;
plot(t, costInfFullDecent.*par.N, t, costContFullDecent.*par.N)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, costInfu, '--')
legend('inf cost', 'cont cost', 'inf cost_u', 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

sgtitle('decentralized control')




figure(3)
plot(t, costInfu, t, costInfFullDecent+costContFullDecent, t, costInfFull+costContFull, t, Celim_per_unit_time*t )
grid on
legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
xlabel('time (days)')
ylabel('cost (units of k)')

