function [results] = getResultsAnalytic(probType, a0, par)

% Iteration settings
maxTries = 1000;
relTol = 1e-8;

% Relaxation factor for iteration - decreasing with iteration count
relFact = 0.3 + 0.7*exp(-15*linspace(0, 1, maxTries));

% Vector of times
t = 0:par.dt:par.tMax;

% Solving problem for centralised or decentralised control?
if probType == "cent"
    coeff = 2;
elseif probType == "decent"
    coeff = 1;
else
    error('probType must be "cent" or "decent')
end

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

    % Calculate analytical solution
    aNew = (par.costlin + 2*par.costquad)./(2*par.costquad + coeff*par.costPerInf.*  (par.Beta*I)./par.N  .*  S(:, end)./par.N  );

    % Update a using specified relaxation factor
    a = max(0, min(1, (1-relFact(iTry))*a + relFact(iTry)*aNew ));

    % Test for convergence
    dx =  max(max(abs(a-aSav))) / max(max(abs(aSav)));
    convFlag = dx <= relTol;
    iTry = iTry+1;
end

if ~convFlag
    fprintf('Waring: analytic method not converged, dx = %.4e\n', dx)

    figure(1);
    hold on

    for ii = 1:4

    
        % Solve SIR model
        [S, I, ~] = solveModelFull(a, a, par);
    
        % Calculate analytical solution
        aNew = (par.costlin + 2*par.costquad)./(2*par.costquad + coeff*par.costPerInf.*  (par.Beta*I)./par.N  .*  S(:, end)./par.N  );
    
        % Update a using specified relaxation factor
        a = max(0, min(1, (1-relFact(end))*a + relFact(end)*aNew ));

       plot(t, a)
       pause
    end
    0
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

