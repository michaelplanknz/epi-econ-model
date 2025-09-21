function [f, costInf, costCont, a, S, I, S_focal] = objFnDecentFull(x, xBG, focalGroupID, par)


t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

a = pchip(tMesh, x', t);
aBG = pchip(tMesh, xBG', t);

% x is just the control parameters for the focal group
% solveModel requires focal and background parameters for both groups so
% just assume the non-focal group is doing the same as the focal group (this should not affect
% output, which only appleis to the focal group)
atmp = repmat(a, par.nGroups, 1);

[S, I, S_focal] = solveModelFull(atmp, aBG, par);

S_focal = S_focal(focalGroupID, :);

costInf = par.costPerInf(focalGroupID)*(1-S_focal/par.N(focalGroupID)); 
costCont = par.dt * cumsum(par.costlin(focalGroupID).*(1-a) + par.costquad(focalGroupID).*(1-a).^2, 2);

f = costInf(end)+costCont(end);



