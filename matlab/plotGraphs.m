clear 
close all

% Set output folder location
outFolder = '../output/';

fIn = outFolder + "results.mat";
load(fIn);


% Main results use cost per infection = 1 and control costs are relative to this
% Use this to scale to all costs to $
dollarsPerInf = 10000;

% Compute elimination costs
[Celim, aElimOut] = calcElimCost(par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = figure(1);
h.Position = [     240         125        1024         760];
tiledlayout(2, 2, "TileSpacing", "compact");

% Centralised problem
nexttile;
plot(t, resultsCent.R/sum(par.N))
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.R/sum(par.N), '--')
plot(t, resultsCent.a, 'LineWidth', 2)
if par.nGroups > 1
    lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
else
    lbls = ["R(t)", "R(t) (unmitigated)", "a(t)"];
end
legend(lbls, 'Location', 'southeast')
grid on
xlabel('time (days)')
title('(a) epidemic dynamics - centralized control')

nexttile;
plot(t, (resultsCent.costInf + resultsCent.costCont)*dollarsPerInf/1e9)
hold on
plot(t, resultsCent.costInf*dollarsPerInf/1e9)
plot(t, resultsCent.costCont*dollarsPerInf/1e9)
set(gca, 'ColorOrderIndex', 1);
plot(t, (results_u.costInf + results_u.costCont)*dollarsPerInf/1e9, '--')
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["total cost", "infection cost", "control cost", "total cost (unmitigated)"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative cost ($bn)')
xlabel('time (days)')
title('(b) costs - centralized control')


% Decentralised problem
nexttile;
plot(t, resultsDecent.R/sum(par.N))
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.R/sum(par.N), '--')
plot(t, resultsDecent.a, 'LineWidth', 2)
grid on
if par.nGroups > 1
    lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
else
    lbls = ["R(t)", "R(t) (unmitigated)", "a(t)"];
end
legend(lbls, 'Location', 'southeast')
xlabel('time (days)')
title('(c) epidemic dynamics - decentralized control')


nexttile;
plot(t, (resultsDecent.costInf + resultsDecent.costCont)*dollarsPerInf/1e9)
hold on
plot(t, resultsDecent.costInf*dollarsPerInf/1e9)
plot(t, resultsDecent.costCont*dollarsPerInf/1e9)
set(gca, 'ColorOrderIndex', 1);
plot(t, (results_u.costInf + results_u.costCont)*dollarsPerInf/1e9, '--')
grid on
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["total cost", "infection cost", "control cost", "total cost (unmitigated)"];
end
legend(lbls, 'Location', 'southeast')
ylabel('cumulative cost ($bn)')
xlabel('time (days)')
title('(d) costs - decentralized control')



% For the purposes of plotting, simulate the elimination response to sporadic
% outbreaks
aElim = ones(size(t));

% Assume outbreaks are evenly spaced with the first one at 0.5 * expected
% spacing
tOut1 = 0.5/par.r;

% Mean outbreak duration
tDur = log(par.xOutbreak)/(par.Gamma-par.Beta*par.alpha_TTI*aElimOut^2);

tRel = mod(t-tOut1, 1/par.r);
aElim(tRel >= 0 & tRel < tDur) = aElimOut;



% Cost comparison
h = figure(3);
h.Position = [703         452        1024         380];
tiledlayout(1, 2, "TileSpacing", "compact")

nexttile;
plot(t, results_u.a)
hold on
plot(t, resultsDecent.a)
plot(t, resultsCent.a)
plot(t, aElim)
ylim([0, 1])
grid on
xlabel('time (days)')
ylabel('a(t)')
legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
title('(a) relative contact rate')

nexttile;
plot(t, results_u.costInf*dollarsPerInf/1e9)
hold on
plot(t, (resultsDecent.costInf + resultsDecent.costCont)*dollarsPerInf/1e9)
plot(t, (resultsCent.costInf + resultsCent.costCont)*dollarsPerInf/1e9)
plot(t, Celim*t*dollarsPerInf/1e9 )
grid on
xlabel('time (days)')
ylabel('cumulative cost ($bn)')
title('(b) costs')

