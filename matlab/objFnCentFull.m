function [f, results] = objFnCentFull(x, par)

t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

% Interpolate the control function a(t) between meshpoints as a piecewise
% cubic
a = max(0, min(1, pchip(tMesh, x, t)));

% Solve SIR model
[S, I, ~] = solveModelFull(a, a, par);

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
results.x = x;

% Calculate objective function
f = norm(costInfPP(:, end)+costContPP(:, end), par.normP);

% Opptional - if under HIT add cost of infecting enough people to reach HIT
% R0 = par.Beta/par.Gamma;
% f = norm(costInfPP(:, end)+costContPP(:, end), par.normP) + par.costPerInf*max(0, S(end)/par.N-1/R0 );


