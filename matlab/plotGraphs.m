clear 
close all

% Set output folder location
outFolder = '../output/';
figFolder = '../figures/';

% Set to true to save Figures as .png files
saveFlag = false;


% Load previously saved model results
fIn = outFolder + "results.mat";
load(fIn);

% Select values of Beta and costPerInf to plot
Beta_arr       = [0.3   0.3   0.3   0.6   0.6   0.6];
costPerInf_arr = [0.1   0.5   1.2   0.1   0.5   1.2];
nPlots = length(Beta_arr);

% Main results use cost per infection measured in $10,000s
% This scaling factor converts all costs to $
dollarsPerInf = 10000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - individual scenario plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iPlot = 1:nPlots

    iScenario = find(Beta_list == Beta_arr(iPlot) & costPerInf_list == costPerInf_arr(iPlot));

    %Set scneario-dependent parameters
    par.Beta = Beta_arr(iPlot);
    par.costPerInf = costPerInf_arr(iPlot);

    % Compute elimination costs
    [Celim, aElimOut] = calcElimCost(par);

    h = figure(2*iPlot-1);
    h.Position = [     240         125        1024         760];
    tiledlayout(2, 2, "TileSpacing", "compact");
    
    % Centralised problem
    nexttile;
    plot(t, resultsCent(iScenario).R/sum(par.N))
    hold on
    set(gca, 'ColorOrderIndex', 1);
    plot(t, results_u(iScenario).R/sum(par.N), '--')
    plot(t, resultsCent(iScenario).a, 'LineWidth', 2)
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
    plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
    hold on
    plot(t, resultsCent(iScenario).costInf*dollarsPerInf/1e9)
    plot(t, resultsCent(iScenario).costCont*dollarsPerInf/1e9)
    set(gca, 'ColorOrderIndex', 1);
    plot(t, (results_u(iScenario).costInf + results_u(iScenario).costCont)*dollarsPerInf/1e9, '--')
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
    sgtitle("R_0=" + par.Beta/par.Gamma + ", k=$" + par.costPerInf*dollarsPerInf)
    
    
    % Decentralised problem
    nexttile;
    plot(t, resultsDecent(iScenario).R/sum(par.N))
    hold on
    set(gca, 'ColorOrderIndex', 1);
    plot(t, results_u(iScenario).R/sum(par.N), '--')
    plot(t, resultsDecent(iScenario).a, 'LineWidth', 2)
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
    plot(t, (resultsDecent(iScenario).costInf + resultsDecent(iScenario).costCont)*dollarsPerInf/1e9)
    hold on
    plot(t, resultsDecent(iScenario).costInf*dollarsPerInf/1e9)
    plot(t, resultsDecent(iScenario).costCont*dollarsPerInf/1e9)
    set(gca, 'ColorOrderIndex', 1);
    plot(t, (results_u(iScenario).costInf + results_u(iScenario).costCont)*dollarsPerInf/1e9, '--')
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
    h = figure(2*iPlot);
    h.Position = [703         452        1024         380];
    tiledlayout(1, 2, "TileSpacing", "compact")
    
    nexttile;
    plot(t, results_u(iScenario).a)
    hold on
    plot(t, resultsDecent(iScenario).a)
    plot(t, resultsCent(iScenario).a)
    plot(t, aElim)
    ylim([0, 1])
    grid on
    xlabel('time (days)')
    ylabel('a(t)')
    legend('unmitigated', 'decentralised', 'centralised', 'elimination', 'Location', 'southeast')
    title('(a) relative contact rate')
    
    nexttile;
    plot(t, results_u(iScenario).costInf*dollarsPerInf/1e9)
    hold on
    plot(t, (resultsDecent(iScenario).costInf + resultsDecent(iScenario).costCont)*dollarsPerInf/1e9)
    plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
    plot(t, Celim*t*dollarsPerInf/1e9 )
    grid on
    xlabel('time (days)')
    ylabel('cumulative cost ($bn)')
    title('(b) costs')
    sgtitle("R_0=" + par.Beta/par.Gamma + ", k=$" + par.costPerInf*dollarsPerInf)

    if saveFlag
        fName = "fig"+iPlot+".png";
        saveas(h, figFolder+fName);
    end

end


if ~exist('Beta_mat')
    [Beta_mat, costPerInf_mat] = meshgrid(Beta_vals, costPerInf_vals);
end

% Initialise matrices for heat map results
costUnmit = zeros(size(Beta_mat));
costDecent = zeros(size(Beta_mat));
costCent = zeros(size(Beta_mat));
costElim = zeros(size(Beta_mat));
HIT_shortfall = zeros(size(Beta_mat));
stratCode = nan(size(Beta_mat));
tCrit = nan(size(Beta_mat));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - heat maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iScenario = 1:nScenarios

    %Set scneario-dependent parameters
    par.Beta = Beta_list(iScenario);
    par.costPerInf = costPerInf_list(iScenario);

    % Get matrix indices for this iScenario
    [iRow, jCol] = ind2sub(size(Beta_mat), iScenario);

    % Fill matrix entries
    costUnmit(iRow, jCol) = results_u(iScenario).costInf(end);
    costDecent(iRow, jCol) = resultsDecent(iScenario).costInf(end) + resultsDecent(iScenario).costInf(end);
    costCent(iRow, jCol) = resultsCent(iScenario).costInf(end) + resultsCent(iScenario).costInf(end);
    HIT_shortfall(iRow, jCol) = resultsCent(iScenario).S(end)/par.N - 1/(par.Beta/par.Gamma);

     % Compute elimination costs
    Celim = calcElimCost(par);
    costElim(iRow, jCol) = Celim*max(t);

    % Calculate threshold time for elimination
    % Note if the centralised response has suppressed the epidemic -
    % threshold hold time is at least max(t) but exact value is unknown -
    % would need to run model longer
    if HIT_shortfall(iRow, jCol) < 0.1
        tCrit(iRow, jCol) = costCent(iRow, jCol)/Celim;
    end

    if costElim(iRow, jCol) < costCent(iRow, jCol)
        stratCode(iRow, jCol) = 3;
    elseif HIT_shortfall(iRow, jCol) > 0.1
        stratCode(iRow, jCol) = 2;
    else
        stratCode(iRow, jCol) = 1;
    end

end

cMax = max(max(costDecent))*dollarsPerInf/1e9;

h = figure(100);
tiledlayout(2, 2, "TileSpacing", "compact");
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costDecent*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('decentralised cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costCent*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('centralised cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costElim*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('elimination cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, stratCode);
%contourf(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, stratCode, [0.5 1.5, 2.5]);
h = gca; h.YDir = 'normal';
h.Colormap = parula;
xlabel('R_0')
ylabel('cost per infection')
title('optimal strategy')

figure(101);
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, tCrit);
colorbar;
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('threshold time (days)')


