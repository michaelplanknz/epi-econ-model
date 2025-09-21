function [f, costInf, costCont, a, S, I, S_focal] = objFnDecent(x, xBG, focalGroupID, par)



% x is just the control parameters for the focal group
% solveModel requires focal and background parameters for both groups so
% just assume the non-focal group is doing the same as the focal group (this should not affect
% output, which only appleis to the focal group)
xtmp = repmat(x, par.nGroups, 1);

[S, I, S_focal] = solveModel(xtmp, xBG, par);

a = myContModel(xtmp, S, I, par);


S_focal = S_focal(focalGroupID, :);
a = a(focalGroupID, :);

costInf = par.costPerInf(focalGroupID)*(1-S_focal/par.N(focalGroupID)); 
costCont = par.dt * cumsum(par.costlin(focalGroupID).*(1-a) + par.costquad(focalGroupID).*(1-a).^2, 2);

f = costInf(end)+costCont(end);



