function [f, results] = objFnCentFull(x, par)

t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

% Interpolate the control function a(t) between meshpoints as a piecewise
% cubic
a = pchip(tMesh, x, t);

% Solve SIR model
[S, I, ~] = solveModelFull(a, a, par);

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



