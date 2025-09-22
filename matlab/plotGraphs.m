clear 
close all

% Set output folder location
outFolder = '../output/';

fIn = outFolder + "results.mat";
load(fIn);


% Plotting
% Centralised problem
h = figure(1);
h.Position = [     240         125        1024         744];
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
plot(t, resultsCentFull.R)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.R, '--')
plot(t, resultsCentFull.a, 'LineWidth', 2)
if par.nGroups > 1
    lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
else
    lbls = ["R", "R_u", "a"];
end
legend(lbls, 'Location', 'southeast')
grid on
xlabel('time (days)')

nexttile;
plot(t, resultsCentFull.costInf + resultsCentFull.costCont)
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
plot(t, resultsCentFull.costInf, t, resultsCentFull.costCont)
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
plot(t, resultsDecentFull.R)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.R, '--')
plot(t, resultsDecentFull.a, 'LineWidth', 2)
grid on
legend('R', 'R_u', 'a', 'Location', 'southeast')
xlabel('time (days)')


nexttile;
plot(t, (resultsDecentFull.costInf + resultsDecentFull.costCont).*par.N)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf + results_u.costCont, '--')
grid on
legend('cost', 'cost_u', 'Location', 'southeast')
ylabel('cumulative aggregate cost')
xlabel('time (days)')


nexttile;
plot(t, resultsDecentFull.costInf.*par.N, t, resultsDecentFull.costCont.*par.N)
hold on
set(gca, 'ColorOrderIndex', 1);
plot(t, results_u.costInf, '--')
legend('inf cost', 'cont cost', 'inf cost_u', 'Location', 'southeast')
grid on
ylabel('cumulative aggregate cost')
xlabel('time (days)')

sgtitle('decentralized control')



% Cost comparison
figure(3)
plot(t, results_u.costInf, t, resultsDecentFull.costInf + resultsDecentFull.costCont, t, resultsCentFull.costInf + resultsCentFull.costCont, t, Celim_per_unit_time*t )
grid on
legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
xlabel('time (days)')
ylabel('cost (units of k)')
