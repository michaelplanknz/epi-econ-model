function [f, results] = objFnDecentFull(x, xBG, focalGroupID, par)


t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

% Interpolate the control function a(t) between meshpoints as a piecewise
% cubic
a = pchip(tMesh, x', t);
aBG = pchip(tMesh, xBG', t);

% x is just the control parameters for the focal group
% solveModel requires focal and background parameters for both groups so
% just assume the non-focal group is doing the same as the focal group (this should not affect
% output, which only appleis to the focal group)
atmp = repmat(a, par.nGroups, 1);

% Solve SIR model
[S, I, S_focal_all] = solveModelFull(atmp, aBG, par);

% Extract S_focal for the specified focal group
S_focal = S_focal_all(focalGroupID, :);

% Calculate infection and control costs
costInf = par.costPerInf(focalGroupID)*(1-S_focal/par.N(focalGroupID)); 
costCont = par.dt * cumsum(par.costlin(focalGroupID).*(1-a) + par.costquad(focalGroupID).*(1-a).^2, 2);

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



