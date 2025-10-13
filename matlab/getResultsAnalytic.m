function results = getResultsAnalytic(probType, par)

maxTries = 100;
relTol = 1e-8;

t = 0:par.dt:par.tMax;

if probType == "cent"
    coeff = 2;
elseif probType == "decent"
    coeff = 1;
else
    error('probType must be "cent" or "decent')
end

% Initialise at uncontrolled state
a = ones(par.nGroups, length(t));

convFlag = false;
iTry = 1;
while iTry <= maxTries * ~convFlag
    aSav = a;

    % Solve SIR model
    [S, I, ~] = solveModelFull(a, a, par);

    % Update a
    a = max(0, min(1, (par.costlin + 2*par.costquad)./(2*par.costquad + coeff*par.costPerInf.*  (par.Beta*I)./par.N  .*  S(:, end)) ));

    % Test for convergence
    convFlag = max(max(abs(a-aSav))) / max(max(abs(aSav))) <= relTol;
    iTry = iTry+1;
end

% Calculate infection and control costs
costInf = par.costPerInf.*(par.N-S);
costCont = par.dt * par.N .* cumsum(par.costlin.*(1-a) + par.costquad.*(1-a).^2, 2);

% Store results in structure for output
results.S = S;
results.I = I;
results.R = par.N-S-I;
results.a = a;
results.costInf = costInf;
results.costCont = costCont;


