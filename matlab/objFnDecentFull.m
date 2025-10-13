function [f, results] = objFnDecentFull(a, aBG, par)

% Solve SIR model
[S, I, S_focal] = solveModelFull(a, aBG, par);


% Calculate infection and control costs per person
costInfPP = par.costPerInf .* (1-S_focal./par.N); 
costContPP = par.dt * cumsum(par.costlin.*(1-a) + par.costquad*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.S_focal = S_focal;
results.a = a;
results.costInf = costInfPP.*par.N;
results.costCont = costContPP.*par.N;

% Calculate objective function
f = norm(costInfPP(:, end)+costContPP(:, end), par.normP);



