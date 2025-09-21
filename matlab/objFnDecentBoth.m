function [f, costInf, costCont, a, S, I, S_focal] = objFnDecentBoth(x, xBG, par)

[S, I, S_focal] = solveModel(x, xBG, par);

a = myContModel(x, S, I, par);

costInf = par.costPerInf.*(1-S_focal./par.N);
costCont = cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

f = costInf(:, end) + costCont(:, end);



