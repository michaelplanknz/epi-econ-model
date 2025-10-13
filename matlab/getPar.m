function par = getPar()

% Number of days to simulate
par.tMax = 400;

par.dt = 0.1;

% Spacing of time mesh points for defining control function (days)
par.meshSpace = 10;



% Recovery rate (days-1)
par.Gamma = 0.2;

% Initial infected fraction
par.I0 = 1e-5;

% Norm power for calculating aggregated cost across mutlipe groups (set to
% 1 to just sum costs across groups)
par.normP = 1;

% Population size (in each group)
par.N = 5e6;
par.nGroups = length(par.N);

% Transmission rate matrix
%par.Beta = 0.3*[1.6 0.4; 0.4 0.1];
par.Beta = 0.4;

% Linear and quadratic cost coefficients in each group
par.costlin = 0;
par.costquad = 0.04;

% Cost per infection in each group
par.costPerInf = 1;



% Elimination model parameters 

% Outbreak size at first detection
par.xOutbreak = 20;

% Outbreak frequency (days-1)
par.r = 1/150;

% Border closure cost per unit time (set to 5% of cost of reducing Reff to
% 1)
R0 = eigs(par.Beta, [], 1)/par.Gamma;
costlin_avg = sum(par.costlin.*par.N);
costquad_avg = sum(par.costquad.*par.N);
par.b = 0.05 * (costlin_avg*(1-sqrt(1/R0)) + costquad_avg*(1-sqrt(1/R0))^2);

% Fraction of country under control measures to eliminate border related
% outbreaks
par.controlFrac = 0.37;

% Multiplicative effect of TTI measures on R0
par.alpha_TTI = 0.75;

