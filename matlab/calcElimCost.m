function [CElim, aOptElim, CSup, aSup] = calcElimCost(par)

% Function to calculatet the elimination cost per unit time according to
% the model

% Calculate linear and quadratic costs per unit time if same control is
% applied to a proportion (controlFrac) of the pop
costlin_avg = par.controlFrac * sum(par.costlin.*par.N);
costquad_avg = par.controlFrac * sum(par.costquad.*par.N);

% Calculate optimal a for elimination of small outbreaks
R0 = eigs(par.Beta, [], 1)/par.Gamma;
myF = @(a)( (costlin_avg*(1-a) + costquad_avg*(1-a).^2)./(par.Gamma - par.Gamma*R0*par.alpha_TTI*a.^2) );

opts = optimset('display', 'off');
aOptElim = fmincon(myF, 0, [], [], [], [], 0, sqrt(1/(par.alpha_TTI*R0)), [], opts );

% Calculate elimination cost per unit time
CElim = par.b + par.r*log(par.xOutbreak)*(costlin_avg*(1-aOptElim) + costquad_avg*(1-aOptElim)^2)/(par.Gamma - R0*par.Gamma*par.alpha_TTI*aOptElim^2);


% Calculate suppression cost per unit time
aSup = sqrt(1/(par.alpha_TTI*R0));
CSup = sum(par.costlin.*par.N)*(1-aSup) + sum(par.costquad.*par.N)*(1-aSup)^2;

