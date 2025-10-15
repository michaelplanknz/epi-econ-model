function a = myContModel(x, S, I, par)

% Heuristic model for relative contact rate (relative to preferred value)
% in terms of FOI and susceptibility S(t)

% Model has two coefficients xa and xb
% If there are multiple goupts these are stored as [xa1 xb1 xa2 xb2 ...]
xa = x(1:2:end);
xb = x(2:2:end);

% Calculate uncontrolled FOI
[aFocal, aBG] = deal( ones(par.nGroups, 1));
FOI = calcFOI(I, aFocal, aBG, par);

% Evaluate heuristic contact function
a = (1-xa.*S./par.N)./(1+xb.*FOI.*S./par.N);

