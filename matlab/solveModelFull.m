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


% 
% 
% nGroups = length(par.N);
% 
% tSpan = 0:par.nDays-1;
% 
% % IC for S, I and S_focal
% S0 = par.N*(1-par.I0);
% I0 = par.N*par.I0;
% S_focal0 = S0;
% 
% IC = [S0; I0; S_focal0];
% opts = odeset('NonNegative', 1:length(IC));
% 
% [t, Y] = ode45(@(t, y)myRHSFull(t, y, x_focal, xBG, par), tSpan, IC);
% 
% S = Y(:, 1:nGroups)';
% I = Y(:, nGroups+1:2*nGroups)';
% S_focal = Y(:, 2*nGroups+1:3*nGroups)';








