function par = getPar()

% Number of days to simulate
par.tMax = 600;

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
par.Beta = 0.3;

% Linear and quadratic cost coefficients in each group (units of $10k)
par.costlin = 0;
par.costquad = 0.02;

% Cost per infection in each group (units of $10k)
par.costPerInf = 0.1;



% Elimination model parameters 

% Outbreak size at first detection
par.xOutbreak = 20;
par.tDet = 14;

% Outbreak frequency (days-1)
par.r = 1/150;

% Border closure cost per unit time
par.b = 3000;

% Fraction of country under control measures to eliminate border related
% outbreaks
par.controlFrac = 0.4;

% Multiplicative effect of TTI measures on R0 (alpha=1 is no effect)
par.alpha_TTI = 0.8;

% Max. prevalence (%) at which TTI works 
par.TTI_max = 1000;

% Breadth of logistic transition of TTI from on to off
par.TTI_breadth = 50;
