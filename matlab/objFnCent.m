function [f, results] = objFnCent(x, par)

% Solve SIR model
[S, I, ~] = solveModel(x, x, par);

% Evaluate heuristic function for contact rate a(t)
a = myContModel(x, S, I, par);

% Calculate infection and control costs
costInf = par.costPerInf.*(par.N-S);
costCont = par.dt * par.N .* cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.a = a;
results.costInf = costInf;
results.costCont = costCont;

% Calculate objective function
f = norm(costInf(:, end)+costCont(:, end), par.normP);



