function [f, costInf, costCont, a, S, I] = objFnCent(x, par)

[S, I, ~] = solveModel(x, x, par);

a = myContModel(x, S, I, par);

costInf = par.costPerInf.*(par.N-S);
costCont = par.dt * par.N .* cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

f = norm(costInf(:, end)+costCont(:, end), par.normP);



