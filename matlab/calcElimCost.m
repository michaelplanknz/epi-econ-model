function [CElimRate, aOptElim, CElimSim, aElimSim, CSupRate, aSup] = calcElimCost(t, par)

% Function to calculatet the elimination and suppression cost per unit time according to
% the model

R0 = eigs(par.Beta, [], 1)/par.Gamma;
dt = t(2)-t(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Suppression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate suppression activity level
aSup = sqrt(1/(par.alpha_TTI*R0));

% Calculate suppression cost per unit time
CSupRate = sum(par.costlin.*par.N)*(1-aSup) + sum(par.costquad.*par.N)*(1-aSup)^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elimination (average)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate linear and quadratic costs per unit time if same control is
% applied to a proportion (controlFrac) of the pop
costlin_avg = par.controlFrac * sum(par.costlin.*par.N);
costquad_avg = par.controlFrac * sum(par.costquad.*par.N);

% Calculate optimal a for elimination of small outbreaks
myF = @(a)( (costlin_avg*(1-a) + costquad_avg*(1-a).^2)./(par.Gamma - par.Gamma*R0*par.alpha_TTI*a.^2) );
opts = optimset('display', 'off');
aOptElim = fmincon(myF, 0, [], [], [], [], 0, aSup, [], opts );

% Calculate elimination cost per unit time
CElimRate = par.b + par.r*par.tDet*(R0-1)*(costlin_avg*(1-aOptElim) + costquad_avg*(1-aOptElim)^2)/(1 - R0*par.alpha_TTI*aOptElim^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Elimination (simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% For the purposes of plotting, simulate the elimination response to sporadic
% outbreaks
aElimSim = ones(size(t));

% Mean outbreak duration
R0 = eigs(par.Beta, [], 1)/par.Gamma;
tDur = par.tDet*(R0-1)/(1 - R0*par.alpha_TTI*aOptElim^2);

% Assume outbreaks are evenly spaced at the expected spacing,
% create a repeating time variable with this frequency
tRel = mod(t, 1/par.r);

% Set aElim to be the elimination level during outbreaks
aElimSim(tRel >= 1/par.r-tDur) = aOptElim;

% Cost rate during outbreak periods
cOut = costlin_avg*(1-aOptElim) + costquad_avg*(1-aOptElim)^2;

% Time-dependent cost rate
y = par.b*ones(size(t));
y(tRel >= 1/par.r-tDur) = par.b + cOut;

% Cumulative cost
CElimSim = dt*cumsum(y);








