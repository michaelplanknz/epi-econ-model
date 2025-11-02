function [f, results] = objFnCentHeur(x, par)

% Evaluate centralised objective function using the heuristic contact rate function

% Solve SIR model with the specified parameters (x) for the heuristic
% contact functio n
[S, I, ~] = solveModelHeur(x, x, par);

% Evaluate heuristic function for contact rate a(t) with the solution for
% S(t) and I(t)
a = myContModel(x, S, I, par);

HIT = 1 - par.Gamma/par.Beta;

% Calculate infection and control costs per person
costInfPP = par.costPerInf.*(1-S./par.N);
costInfPP(:, end) = par.costPerInf.*max(HIT, 1-S(:, end)./par.N);
costContPP = par.dt * cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.a = a;
results.costInf = costInfPP.*par.N;
results.costCont = costContPP.*par.N;
results.x = x;

% Calculate objective function
f = norm(costInfPP(:, end) + costContPP(:, end), par.normP);



