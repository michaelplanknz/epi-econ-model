clear 
close all

% Set output folder location
outFolder = '../output/';

fIn = outFolder + "results.mat";
load(fIn);



% Compute elimination costs
Celim = calcElimCost(par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Centralised problem
h = figure(1);
h.Position = [     240         125        1024         744];
tiledlayout(2, 2, "TileSpacing", "compact");
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
plot(t, resultsCent.costInf + resultsCent.costCont)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf + results_u.costCont, '--')
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["cost", "cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

nexttile;
plot(t, resultsCent.costInf, t, resultsCent.costCont)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf, '--')
if par.nGroups > 1
    lbls = ["inf cost_"+groupIDs; "cont cost_"+groupIDs; "inf cost_"+groupIDs+"_u"];
else
    lbls = ["inf cost", "cont cost", "inf cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

sgtitle('centralized control')
drawnow


% Decentralised problem
h = figure(2);
h.Position = [    360         183        1024         755];
tiledlayout(2, 2, "TileSpacing", "compact");
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
plot(t, resultsDecent.costInf + resultsDecent.costCont)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf + results_u.costCont, '--')
grid on
if par.nGroups > 1
    lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
else
    lbls = ["cost", "cost_u"];
end
legend(lbls, 'Location', 'southeast')
ylabel('cumulative aggregate cost')
xlabel('time (days)')


nexttile;
plot(t, resultsDecent.costInf, t, resultsDecent.costCont)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf, '--')
if par.nGroups > 1
    lbls = ["inf cost_"+groupIDs; "cont cost_"+groupIDs; "inf cost_"+groupIDs+"_u"];
else
    lbls = ["inf cost", "cont cost", "inf cost_u"];
end
legend(lbls, 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

sgtitle('decentralized control')



% Cost comparison
figure(3)
plot(t, results_u.costInf, t, resultsDecent.costInf + resultsDecent.costCont, t, resultsCent.costInf + resultsCent.costCont, t, Celim*t )
grid on
legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
xlabel('time (days)')
ylabel('cost (units of k)')
