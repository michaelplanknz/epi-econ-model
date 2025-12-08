clear 
close all

% Set output folder location
outFolder = '../output/';
figFolder = '../figures/';

% Set to true to save Figures as .png files
saveFlag = false;


% Select values of Beta and costPerInf to plot
Beta_arr       = [0.3, 0.6];
costPerInf_arr = [0.2   0.6   1.2];
nPlots = length(Beta_arr);
nSubplots = length(costPerInf_arr);

% Load previously saved model results
fIn = outFolder + "results_for_plots.mat";
load(fIn);

Beta_arr       = [0.3, 0.6];
costPerInf_arr = [0.2   0.6   1.2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - individual scenario plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

letters = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)"   ];

for iPlot = 1:nPlots
    h = figure(iPlot);
    h.Position = [139          69        1024         912];
    tiledlayout(nSubplots, 2, "TileSpacing", "compact")

    h = figure(iPlot+nPlots);
    h.Position = [139          69        1024         912];
    tiledlayout(nSubplots, 2, "TileSpacing", "compact")
    
    h = figure(iPlot+2*nPlots);
    h.Position = [139          69        1500         350];
    tiledlayout(1, nSubplots, "TileSpacing", "compact")


    for iSubplot = 1:nSubplots

        iScenario = find(Beta_list == Beta_arr(iPlot) & costPerInf_list == costPerInf_arr(iSubplot));
    
        %Set scneario-dependent parameters
        par.Beta = Beta_arr(iPlot);
        par.costPerInf = costPerInf_arr(iSubplot);
    
        % Compute elimination costs
        [CElimRate, ~, CElimSim, aElimSim, CSupRate, aSup] = calcElimCost(t, par);
    

        % Calculate dC/da as a check
        Sinf = resultsDecent(iScenario).S(end)/par.N;
        dCda = par.costPerInf*par.Beta/par.N * Sinf * resultsDecent(iScenario).a .* resultsDecent(iScenario).I - par.costlin - 2*par.costquad*(1-resultsDecent(iScenario).a);
        indControl = find(resultsDecent(iScenario).a < 0.99999999);

        % Main scenario plot
        figure(iPlot);
        nexttile;
        plot(t, results_u(iScenario).a.^2, 'LineWidth', 2)
        hold on
        plot(t, resultsDecent(iScenario).a.^2, 'LineWidth', 2)
        plot(t, resultsCent(iScenario).a.^2, 'LineWidth', 2)
        plot(t, aSup^2*ones(size(t)), 'LineWidth', 2)
        plot(t, aElimSim.^2, 'LineWidth', 2)
        set(gca, 'ColorOrderIndex', 1);
        plot(t, results_u(iScenario).R/sum(par.N), '--')
        plot(t, resultsDecent(iScenario).R/sum(par.N), '--')
        plot(t, resultsCent(iScenario).R/sum(par.N), '--')
        ylim([0, 1])
        xlim([0 tHoriz])
        grid on
        xlabel('time (days)')
        ylabel('[cumulative infections  ,  a(t)^2]')
        title(letters(2*iSubplot-1) + " cost per infection = $" + dollarsPerInf*par.costPerInf)
 
        % Cost comparison
        nexttile;
        plot(t, results_u(iScenario).costInf*dollarsPerInf/1e9)
        hold on
        plot(t, (resultsDecent(iScenario).costInf + resultsDecent(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, CSupRate*t*dollarsPerInf/1e9 )
        plot(t, CElimSim*dollarsPerInf/1e9 )
%         set(gca, 'ColorOrderIndex', 1);
%         plot(t, CElimRate*t*dollarsPerInf/1e9, ':')
        xlim([0 tHoriz])
        ylim([0 cUpper(iPlot)])
        grid on
        xlabel('time (days)')
        ylabel('cumulative cost ($bn)')
        title(letters(2*iSubplot) + " cost per infection = $" + dollarsPerInf*par.costPerInf)



        % State-dependent version of the same figure
        figure(iPlot+nPlots);
        nexttile;
        plot(t, results_u(iScenario).a.^2, 'LineWidth', 2)
        hold on
        plot(t, resultsDecent_SD(iScenario).a.^2, 'LineWidth', 2)
        plot(t, resultsCent(iScenario).a.^2, 'LineWidth', 2)
        plot(t, aSup^2*ones(size(t)), 'LineWidth', 2)
        plot(t, aElimSim.^2, 'LineWidth', 2)
        set(gca, 'ColorOrderIndex', 1);
        plot(t, results_u(iScenario).R/sum(par.N), '--')
        plot(t, resultsDecent_SD(iScenario).R/sum(par.N), '--')
        plot(t, resultsCent(iScenario).R/sum(par.N), '--')
        ylim([0, 1])
        xlim([0 tHoriz])
        grid on
        xlabel('time (days)')
        ylabel('[cumulative infections  ,  a(t)^2]')
        title(letters(2*iSubplot-1) + " cost per infection = $" + dollarsPerInf*par.costPerInf)
 
        % Cost comparison
        nexttile;
        plot(t, results_u(iScenario).costInf*dollarsPerInf/1e9)
        hold on
        plot(t, (resultsDecent_SD(iScenario).costInf + resultsDecent_SD(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, CSupRate*t*dollarsPerInf/1e9 )
        plot(t, CElimSim*dollarsPerInf/1e9 )
%         set(gca, 'ColorOrderIndex', 1);
%         plot(t, CElimRate*t*dollarsPerInf/1e9, ':')
        xlim([0 tHoriz])
        ylim([0 cUpper(iPlot)])
        grid on
        xlabel('time (days)')
        ylabel('cumulative cost ($bn)')
        title(letters(2*iSubplot) + " cost per infection = $" + dollarsPerInf*par.costPerInf)


        figure(iPlot+2*nPlots);
        nexttile;
        plot(t, dCda*dollarsPerInf)
        hold on
        if length(indControl) >= 2
            fill(t(indControl([1, end, end, 1])), [-50, -50, 0, 0 ], [0.4 0.4 0.4], 'FaceAlpha', 0.2, 'LineStyle', 'none'  )
        end
        xlim([0 tHoriz])
        ylim([-50, 0])
        grid on
        xlabel('time (days)')
        ylabel('\delta C/\delta a ($ per person per day)')
        title(letters(iSubplot) + " cost per infection = $" + dollarsPerInf*par.costPerInf)
    end
    
    h = figure(iPlot);
    sgtitle("R_0=" + par.Beta/par.Gamma)
    l = legend('unmitigated', 'mitigation (decent.)', 'mitigation (cent.)', 'suppression', 'elimination', 'Location', 'southoutside');
    l.Layout.Tile = 'south';
    if saveFlag
        fName = "fig" + iPlot + ".png";
        saveas(h, figFolder + fName);
    end

    
    h = figure(iPlot+nPlots);
    sgtitle("R_0=" + par.Beta/par.Gamma)
    l = legend('unmitigated', 'mitigation (state-dep decent.)', 'mitigation (cent.)', 'suppression', 'elimination', 'Location', 'southoutside');
    l.Layout.Tile = 'south';
    if saveFlag
        fName = "figS" + iPlot + ".png";
        saveas(h, figFolder + fName);
    end


    h = figure(iPlot+2*nPlots);
    sgtitle("R_0=" + par.Beta/par.Gamma)
    if saveFlag
        fName = "figS" + (iPlot+nPlots) + ".png";
        saveas(h, figFolder + fName);
    end 
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - heat maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cMax = max(max(costDecent))*dollarsPerInf/1e9;

% Record the threshold value of R0 above hwich elimination is better than
% suppression
ind1 = find(costElim(1, :) > costSup(1, :), 1, 'last');
ind2 = find(costElim(1, :) <= costSup(1, :), 1, 'first');
R0crit = mean(Beta_vals([ind1 ind2]))/par.Gamma;

iPlot = 3*nPlots + 1;
iFile = nPlots+1;
h = figure(iPlot);
h.Position = [   141   407   770   561];
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costDecent*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection ($)')
title('(a) mitigation (decent.) cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costMit*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection ($)')
title('(b) mitigation (cent.) cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, min(costElim, costSup)*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xline(R0crit, 'w--')
xlabel('R_0')
ylabel('cost per infection ($)')
title('(c) elimination/suppression cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, stratCode);
h = gca; h.YDir = 'normal';
h.Colormap = parula;
xlabel('R_0')
text(1.4, 8000, 'suppression', 'Rotation', 90)
text(2.4, 12000, 'elimination')
text(2.2, 5000, 'mitigation', 'Color', 'w')
text(2.2, 1200, 'no PSHMs', 'Color', 'w')
ylabel('cost per infection ($)')
title('(d) optimal strategy')

if saveFlag
    fName = "fig" + iFile + ".png";
    saveas(h, figFolder + fName);
end



iPlot = iPlot + 1;
iFile = iFile + 1;
h = figure(iPlot);
h.Position = [   139   463   994   385];
tiledlayout(1, 2, "TileSpacing", "compact");

nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, tCrit);
cb = colorbar;
h = gca; h.YDir = 'normal';
h.Colormap = hot;
cb.Label.String = 'threshold time (days)';
cb.Label.Rotation = 270;
cb.Label.Position(1) = 4;
xlabel('R_0')
ylabel('cost per infection ($)')
title('(a)')

nexttile;
imagesc(tDetRet_vals, alpha_TTI_vals, tCrit2);
cb = colorbar;
h = gca; h.YDir = 'normal';
h.Colormap = hot;
cb.Label.String = 'threshold time (days)';
cb.Label.Rotation = 270;
cb.Label.Position(1) = 3.8;
xlabel('border outbreak detect-return time ratio')
ylabel('\alpha_{TTI}')
title('(b)')


if saveFlag
    fName = "fig" + iFile + ".png";
    saveas(h, figFolder + fName);
end



iPlot = iPlot + 1;
iFile = iFile + 1;
h = figure(iPlot);
h.Position = [  216   409   994   385];
tiledlayout(1, 2, 'TileSpacing', "compact")
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, (costDecent-costMit)*dollarsPerInf/1e9);
colorbar;
clim([-1 8]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection ($)')
title('(b) state-independent decentralised model')

nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, (costDecent_SD-costMit)*dollarsPerInf/1e9);
colorbar;
clim([-1 8]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection ($)')
title('(b) state-dependent decentralised model')

sgtitle('difference in cost between decentralised and centralised ($ bn)')


if saveFlag
    fName = "figS5.png";
    saveas(h, figFolder + fName);
end






% Extra figure showing the (ignored) infection cost of elimination (arising
% from small border-related outbreaks)
iPlot = iPlot + 1;
h = figure(iPlot);
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costInfElim*dollarsPerInf/1e9);
colorbar;
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection ($)')
title('ignored infection cost of elimination ($ bn)')

