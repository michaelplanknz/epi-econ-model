function results = getResultsAnalytic_SD(a0, par)

% Iteration settings
maxTries = 2000;
relTol = 1e-8;

% Relaxation factor for iteration - decreasing with iteration count
relFact = 0.1 + 0.9*exp(-15*linspace(0, 1, maxTries));

% Vector of times
t = 0:par.dt:par.tMax;
nSteps = length(t);

% Set initial condition for a(t)
a = a0;

% Start iterating
convFlag = false;
iTry = 1;
while iTry <= maxTries & ~convFlag
    % Store existing a(t) for convergence testing
    aSav = a;

    % Solve SIR model
    [S, I, ~] = solveModelFull(a, a, par);

    Q = 0;
    for iStep = 1:nSteps
        ti = nSteps+1-iStep;
        
        % Calculate analytical solution
        aNew = (par.costlin + 2*par.costquad).*(S(ti)+I(ti)) ./ (2*par.costquad.*(S(ti)+I(ti)) + (par.Beta*I(ti))./par.N .* (par.costPerInf.*S(end) - Q));
        
        % Update a using specified relaxation factor
        a(ti) = max(0, min(1, (1-relFact(iTry))*a(ti) + relFact(iTry)*aNew ));

        Q = Q + (S(ti)+I(ti)) * (par.costlin*(1-a(ti)) + par.costquad*(1-a(ti))^2) * par.dt;
    end

    % Test for convergence
    dx =  max(max(abs(a-aSav))) / max(max(abs(aSav)));
    convFlag = dx <= relTol;
    iTry = iTry+1;
end

if ~convFlag
    fprintf('Waring: SD analytic method not converged, dx = %.4e\n', dx)    
end




% Calculate infection and control costs per person
costInfPP = par.costPerInf.*(1-S/par.N);
costContPP = par.dt .* cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.a = a;
results.costInf = costInfPP .* par.N;
results.costCont = costContPP .* par.N;

