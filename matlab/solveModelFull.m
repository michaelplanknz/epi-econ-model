function [S, I, S_focal] = solveModelFull(a_focal, aBG, par)


t = 0:par.dt:par.tMax;
nSteps = length(t)-1;

[S, I, S_focal] = deal(zeros(par.nGroups, nSteps+1));

S(:, 1) = par.N*(1-par.I0);
I(:, 1) = par.N*par.I0;
S_focal(:, 1) = par.N*(1-par.I0);

for iStep = 1:nSteps

    FOI = calcFOI( I(:, iStep), aBG(:, iStep), aBG(:, iStep), par );
    FOI_focal = calcFOI( I(:, iStep), a_focal(:, iStep), aBG(:, iStep), par  );


    S(:, iStep+1) = S(:, iStep).*exp(-FOI*par.dt);
    I(:, iStep+1) = I(:, iStep)*exp(-par.Gamma*par.dt) + S(:, iStep).*(1-exp(-FOI*par.dt));

    S_focal(:, iStep+1) = S_focal(:, iStep).*exp(-FOI_focal*par.dt);
end




