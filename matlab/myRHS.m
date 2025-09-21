function dydt = myRHS(t, y, x, xBG, par)

nGroups = length(par.N);

S = y(1:nGroups);
I = y(nGroups+1:2*nGroups);
S_focal = y(2*nGroups+1:3*nGroups);


a_focal = myContModel(x, S, I, par);
aBG = myContModel(xBG, S, I, par);

FOI = calcFOI( I, aBG, aBG, par );
FOI_focal = calcFOI( I, a_focal, aBG, par  );

dSdt = -FOI.*S;
dIdt = FOI.*S - par.Gamma*I;
dS_focaldt = -FOI_focal.*S_focal;

dydt = [dSdt; dIdt; dS_focaldt];
