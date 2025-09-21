function dydt = myRHSFull(t, y, x_focal, xBG, par)

nGroups = length(par.N);

tMesh = 0:par.meshSpace:par.nDays-1;

S = y(1:nGroups);
I = y(nGroups+1:2*nGroups);
S_focal = y(2*nGroups+1:3*nGroups);

a_focal = pchip(tMesh, x_focal, t);
aBG = pchip(tMesh, xBG, t);

FOI = calcFOI( I, aBG, aBG, par );
FOI_focal = calcFOI( I, a_focal, aBG, par  );

dSdt = -FOI.*S;
dIdt = FOI.*S - par.Gamma*I;
dS_focaldt = -FOI_focal.*S_focal;

dydt = [dSdt; dIdt; dS_focaldt];
