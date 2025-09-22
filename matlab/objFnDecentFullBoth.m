function [f, results] = objFnDecentFullBoth(x, xBG, par)

% As objFnDecent, but takes control parameters x for both/all groups and returns the model results for both/all groups

t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

% Interpolate the control function a(t) between meshpoints as a piecewise
% cubic
a = pchip(tMesh, x', t);
aBG = pchip(tMesh, xBG', t);


% Solve SIR model
[S, I, S_focal] = solveModelFull(a, aBG, par);


% Calculate infection and control costs
costInf = par.costPerInf.*(1-S_focal./par.N); 
costCont = par.dt * cumsum(par.costlin.*(1-a) + par.costquad*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.S_focal = S_focal;
results.a = a;
results.costInf = costInf;
results.costCont = costCont;

% Calculate objective function
f = costInf(end)+costCont(end);



