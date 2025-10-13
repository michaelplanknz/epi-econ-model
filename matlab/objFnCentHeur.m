function [f, results] = objFnCentHeur(x, par)

% Solve SIR model
[S, I, ~] = solveModel(x, x, par);

% Evaluate heuristic function for contact rate a(t)
a = myContModel(x, S, I, par);

% Calculate infection and control costs per person
costInfPP = par.costPerInf.*(1-S./par.N);
costContPP = par.dt * cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.a = a;
results.costInf = costInfPP.*par.N;
results.costCont = costContPP.*par.N;

% Calculate objective function
f = norm(costInfPP(:, end)+costContPP(:, end), par.normP);



