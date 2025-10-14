function [results] = getResultsAnalytic(probType, a0, par)

maxTries = 100;
relTol = 1e-10;

t = 0:par.dt:par.tMax;

if probType == "cent"
    coeff = 2;
elseif probType == "decent"
    coeff = 1;
else
    error('probType must be "cent" or "decent')
end

a = a0;

convFlag = false;
iTry = 1;
relFact = 0.5;
while iTry <= maxTries & ~convFlag
    aSav = a;

    % Solve SIR model
    [S, I, ~] = solveModelFull(a, a, par);

    % Calculate analytical solution
    aNew = (par.costlin + 2*par.costquad)./(2*par.costquad + coeff*par.costPerInf.*  (par.Beta*I)./par.N  .*  S(:, end)./par.N  );

    % Update a
    a = max(0, min(1, (1-relFact)*a + relFact*aNew ));

    % Test for convergence
    convFlag = max(max(abs(a-aSav))) / max(max(abs(aSav))) <= relTol;
    iTry = iTry+1;
end

if ~convFlag
    fprintf('Waring: analytic method not converged\n')
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


% 
% dC1da_an = par.costPerInf*par.Beta * a.* I/par.N .* S(end)/par.N*par.dt;
% dC2da_an = (- par.costlin - 2*par.costquad*(1-a))  *par.dt;
% 
% C1 = costInfPP(end);
% C2 = costContPP(end);
% da = -0.0001;
% 
% 
% 
% for ia = 1:length(a)
%     aPert = a;
%     aPert(ia) = aPert(ia)+da;
%     [Sp, Ip, Sp_focal] = solveModelFull(aPert, a, par);
%     CPert1 = par.costPerInf*(1-Sp_focal(end)/par.N);
%     CPert2 = par.dt*sum(par.costlin.*(1-aPert) + par.costquad*(1-aPert).^2, 2);
%     dCd1a_num(ia) = (CPert1-C1)/da;
%     dCd2a_num(ia) = (CPert2-C2)/da;
% end
% 
% 
% figure;
% plot(t, dC1da_an, t, dCd1a_num)
% 
% figure;
% plot(t, dC2da_an, t, dCd2a_num)