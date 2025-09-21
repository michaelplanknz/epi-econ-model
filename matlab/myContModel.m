function a = myContModel(x, S, I, par)

% Parametric model for relative contact rate (relative to preferred value)
% in terms of infection prevlanece I(t)


xa = x(1:2:end);
xb = x(2:2:end);

% Unctonrolled FOI
[aFocal, aBG] = deal( ones(par.nGroups, 1));
FOI = calcFOI(I, aFocal, aBG, par);

a = (1-xa.*S)./(1+xb.*FOI.*S);

