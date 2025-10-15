clear 
close all

% Set output folder location
outFolder = '../output/';
figFolder = '../figures/';

% Set to true to save Figures as .png files
saveFlag = true;


% Load previously saved model results
fIn = outFolder + "results.mat";
load(fIn);


% Main results use cost per infection measured in $10,000s
% This scaling factor converts all costs to $
dollarsPerInf = 10000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iScenario = 1:nScenarios

    %Set scneario-dependent parameters
    par.Beta = Beta_arr(iScenario);
    par.costPerInf = costPerInf_arr(iScenario);

    % Compute elimination costs
    [Celim, aElimOut] = calcElimCost(par);

    h = figure(2*iScenario-1);
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
    h = figure(2*iScenario);
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
    legend('unmitigated', 'decentralised response', 'centralised response', 'elimination response', 'Location', 'southeast')
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
    
    if saveFlag
        fName = "fig"+iScenario+".png";
        saveas(h, figFolder+fName);
    end

end