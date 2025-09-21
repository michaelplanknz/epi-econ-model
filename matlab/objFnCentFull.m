function [f, costInf, costCont, a, S, I] = objFnCentFull(x, par)

t = 0:par.dt:par.tMax;
tMesh = 0:par.meshSpace:par.tMax;

a = pchip(tMesh, x, t);

[S, I, ~] = solveModelFull(a, a, par);

costInf = par.costPerInf.*(par.N-S);
costCont = par.dt * par.N .* cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

f = norm(costInf(:, end)+costCont(:, end), par.normP);



