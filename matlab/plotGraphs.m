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

% Centralised problem
h = figure(1);
h.Position = [     240         125        1024         380];
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
plot(t, resultsCent.R/sum(par.N))
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.R/sum(par.N), '--')
plot(t, resultsCent.a, 'LineWidth', 2)
if par.nGroups > 1
    lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
else
    lbls = ["R", "R_u", "a"];
end
legend(lbls, 'Location', 'southeast')
grid on
xlabel('time (days)')

nexttile;
plot(t, (resultsCent.costInf + resultsCent.costCont)*dollarsPerInf/1e9)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, (results_u.costInf + results_u.costCont)*dollarsPerInf/1e9, '--')
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["cost", "cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost ($b)')
xlabel('time (days)')

% nexttile;
% plot(t, resultsCent.costInf, t, resultsCent.costCont)
% hold on
% set(gca, 'ColorOrderIndex', 1);
% plot(t, results_u.costInf, '--')
% if par.nGroups > 1
%     lbls = ["inf cost_"+groupIDs; "cont cost_"+groupIDs; "inf cost_"+groupIDs+"_u"];
% else
%     lbls = ["inf cost", "cont cost", "inf cost_u"];
% end
% legend(lbls, 'Location', 'southeast')
% grid on
% ylabel('cumulative aggregate cost')
% xlabel('time (days)')

sgtitle('centralized control')
drawnow


% Decentralised problem
h = figure(2);
h.Position = [    360         183        1024         380];
tiledlayout(1, 2, "TileSpacing", "compact");
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
    lbls = ["R", "R_u", "a"];
end
legend(lbls, 'Location', 'southeast')
xlabel('time (days)')


nexttile;
plot(t, (resultsDecent.costInf + resultsDecent.costCont)*dollarsPerInf/1e9)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, (results_u.costInf + results_u.costCont)*dollarsPerInf/1e9, '--')
grid on
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["cost", "cost_u"];
end
legend(lbls, 'Location', 'southeast')
ylabel('cumulative aggregate cost ($b)')
xlabel('time (days)')


% nexttile;
% plot(t, resultsDecent.costInf, t, resultsDecent.costCont)
% hold on
% set(gca, 'ColorOrderIndex', 1);
% plot(t, results_u.costInf, '--')
% if par.nGroups > 1
%     lbls = ["inf cost_"+groupIDs; "cont cost_"+groupIDs; "inf cost_"+groupIDs+"_u"];
% else
%     lbls = ["inf cost", "cont cost", "inf cost_u"];
% end
% legend(lbls, 'Location', 'southeast')
% grid on
% ylabel('cumulative aggregate cost')
% xlabel('time (days)')

sgtitle('decentralized control')



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
plot(t, results_u.costInf*dollarsPerInf/1e9)
hold on
plot(t, (resultsDecent.costInf + resultsDecent.costCont)*dollarsPerInf/1e9)
plot(t, (resultsCent.costInf + resultsCent.costCont)*dollarsPerInf/1e9)
plot(t, Celim*t*dollarsPerInf/1e9 )
grid on
xlabel('time (days)')
ylabel('cost ($b)')
legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
nexttile
plot(t, results_u.a)
hold on
plot(t, resultsDecent.a)
plot(t, resultsCent.a)
plot(t, aElim)
ylim([0, 1])
grid on
xlabel('time (days)')
ylabel('a(t)')
