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
Beta_arr       = [0.3, 0.6];
costPerInf_arr = [0.1   0.5   1.2];
nPlots = length(Beta_arr);
nSubplots = length(costPerInf_arr);

% Upper limit for y axis for cost plots for each beta
cUpper = [40 80];

% Main results use cost per infection measured in $10,000s
% This scaling factor converts all costs to $
dollarsPerInf = 10000;


% Time horizon to use for coast comparisons and plots (need to be <= par.tMax
tHoriz = 600;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting - individual scenario plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

letters = [ "(a)", "(b)", "(c)", "(d)", "(e)", "(f)"   ];

for iPlot = 1:nPlots
    h = figure(2*iPlot);
    h.Position = [139          69        1024         912];
    tiledlayout(nSubplots, 2, "TileSpacing", "compact")
     
    for iSubplot = 1:nSubplots

        iScenario = find(Beta_list == Beta_arr(iPlot) & costPerInf_list == costPerInf_arr(iSubplot));
    
        %Set scneario-dependent parameters
        par.Beta = Beta_arr(iPlot);
        par.costPerInf = costPerInf_arr(iPlot);
    
        % Compute elimination costs
        [CElim, aElimOut, CSup, aSup] = calcElimCost(par);

    
        
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
        
   
        
        
        nexttile;
        plot(t, results_u(iScenario).a.^2)
        hold on
        plot(t, resultsDecent(iScenario).a.^2)
        plot(t, resultsCent(iScenario).a.^2)
        plot(t, aSup^2*ones(size(t)))
        plot(t, aElim.^2)
        set(gca, 'ColorOrderIndex', 1);
        plot(t, results_u(iScenario).R/sum(par.N), '--')
        plot(t, resultsDecent(iScenario).R/sum(par.N), '--')
        plot(t, resultsCent(iScenario).R/sum(par.N), '--')
        ylim([0, 1])
        xlim([0 tHoriz])
        grid on
        xlabel('time (days)')
        ylabel('attack rate  ;  R_{EI}(t)')
        title(letters(2*iPlot-1))
 
        % Cost comparison
        nexttile;
        plot(t, results_u(iScenario).costInf*dollarsPerInf/1e9)
        hold on
        plot(t, (resultsDecent(iScenario).costInf + resultsDecent(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
        plot(t, CSup*t*dollarsPerInf/1e9 )
        plot(t, CElim*t*dollarsPerInf/1e9 )
        xlim([0 tHoriz])
        ylim([0 cUpper(iPlot)])
        grid on
        xlabel('time (days)')
        ylabel('cumulative cost ($bn)')
        title(letters(2*iPlot))
    end
    sgtitle("R_0=" + par.Beta/par.Gamma)
    l = legend('unmitigated', 'decentralised', 'mitigation', 'suppression', 'elimination', 'Location', 'southoutside');
    l.Layout.Tile = 'south';


    if saveFlag
        fName = "fig" + iPlot + ".png";
        saveas(h, figFolder + fName);
    end

end



% Initialise matrices for heat map results
costUnmit = zeros(size(Beta_mat));
costDecent = zeros(size(Beta_mat));
costMit = zeros(size(Beta_mat));
costSup = zeros(size(Beta_mat));
costElim = zeros(size(Beta_mat));
%HIT_shortfall = zeros(size(Beta_mat));
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
    costDecent(iRow, jCol) = resultsDecent(iScenario).costInf(end) + resultsDecent(iScenario).costCont(end);
    costMit(iRow, jCol) = resultsCent(iScenario).costInf(end) + resultsCent(iScenario).costCont(end);
    %HIT_shortfall(iRow, jCol) = resultsCent(iScenario).S(end)/par.N - 1/(par.Beta/par.Gamma);

     % Compute elimination costs
    [CElim, ~, CSup] = calcElimCost(par);
    costElim(iRow, jCol) = CElim*tHoriz;
    costSup(iRow, jCol) = CSup*tHoriz;

    % Record code for optimal strategy: 1 = mitigation, 2 = suppression, 3 = elimination
    [~, stratCode(iRow, jCol)] = min([costMit(iRow, jCol), costSup(iRow, jCol), costElim(iRow, jCol)]);


    % Calculate threshold time for elimination
    % Note if the centralised response has suppressed the epidemic -
    % threshold hold time is at least max(t) but exact value is unknown -
    % would need to run model longer
%    if HIT_shortfall(iRow, jCol) < 0.1
        tCrit(iRow, jCol) = costMit(iRow, jCol)/min(CElim, CSup);
%    end


end

cMax = max(max(costDecent))*dollarsPerInf/1e9;

h = figure(100);
h.Position = [   141   407   770   561];
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
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, costMit*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('migitgation cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, min(costElim, costSup)*dollarsPerInf/1e9);
colorbar;
clim([0 cMax]);
h = gca; h.YDir = 'normal';
h.Colormap = hot;
xlabel('R_0')
ylabel('cost per infection')
title('elimination/suppression cost ($ bn)')
nexttile;
imagesc(Beta_vals/par.Gamma, costPerInf_vals*dollarsPerInf, stratCode);
h = gca; h.YDir = 'normal';
h.Colormap = parula;
xlabel('R_0')
text(1.25, 12000, 'suppression')
text(2.3, 12000, 'elimination')
text(2.2, 5000, 'mitigation', 'Color', 'w')
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


