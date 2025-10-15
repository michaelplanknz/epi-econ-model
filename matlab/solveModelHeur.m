function [S, I, S_focal] = solveModelHeur(x, xBG, par)

% Solve model with heuristic contact function a(t)
% x specifies the parameters of the contact function for a focal individual
% xBG specifies the parameters of the backgrojund contact function for the rest of the pop

t = 0:par.dt:par.tMax;
nSteps = length(t)-1;

[S, I, S_focal] = deal(zeros(par.nGroups, nSteps+1));

S(:, 1) = par.N*(1-par.I0);
I(:, 1) = par.N*par.I0;
S_focal(:, 1) = par.N*(1-par.I0);

for iStep = 1:nSteps
    % Evaluate the heuristic contact function with the specified
    % parameters (x) and the current values of S(t) and I(t)
    a_focal = myContModel(x, S(:, iStep), I(:, iStep), par);
    aBG = myContModel(xBG, S(:, iStep), I(:, iStep), par);
    
    % Calculate force of infection
    FOI = calcFOI( I(:, iStep), aBG, aBG, par );
    FOI_focal = calcFOI( I(:, iStep), a_focal, aBG, par  );

    % Solve model equations assuming a fixed FOI for the time step dt
    S(:, iStep+1) = S(:, iStep).*exp(-FOI*par.dt);
    I(:, iStep+1) = I(:, iStep)*exp(-par.Gamma*par.dt) + S(:, iStep).*(1-exp(-FOI*par.dt));

    S_focal(:, iStep+1) = S_focal(:, iStep).*exp(-FOI_focal*par.dt);
end






