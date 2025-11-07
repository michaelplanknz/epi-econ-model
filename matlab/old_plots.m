    % h = figure(2*iPlot-1);
    % h.Position = [     240         125        1024         760];
    % tiledlayout(2, 2, "TileSpacing", "compact");
    % 
    % % Centralised problem
    % nexttile;
    % plot(t, resultsCent(iScenario).R/sum(par.N))
    % hold on
    % set(gca, 'ColorOrderIndex', 1);
    % plot(t, results_u(iScenario).R/sum(par.N), '--')
    % plot(t, resultsCent(iScenario).a, 'LineWidth', 2)
    % if par.nGroups > 1
    %     lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
    % else
    %     lbls = ["R(t)", "R(t) (unmitigated)", "a(t)"];
    % end
    % legend(lbls, 'Location', 'southeast')
    % xlim([0 tHoriz])
    % grid on
    % xlabel('time (days)')
    % title('(a) epidemic dynamics - centralized control')
    % 
    % nexttile;
    % plot(t, (resultsCent(iScenario).costInf + resultsCent(iScenario).costCont)*dollarsPerInf/1e9)
    % hold on
    % plot(t, resultsCent(iScenario).costInf*dollarsPerInf/1e9)
    % plot(t, resultsCent(iScenario).costCont*dollarsPerInf/1e9)
    % set(gca, 'ColorOrderIndex', 1);
    % plot(t, (results_u(iScenario).costInf + results_u(iScenario).costCont)*dollarsPerInf/1e9, '--')
    % if par.nGroups > 1
    %     lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
    % else
    %     lbls = ["total cost", "infection cost", "control cost", "total cost (unmitigated)"];
    % end
    % legend(lbls, 'Location', 'southeast')
    % xlim([0 tHoriz])
    % grid on
    % ylabel('cumulative cost ($bn)')
    % xlabel('time (days)')
    % title('(b) costs - centralized control')
    % sgtitle("R_0=" + par.Beta/par.Gamma + ", k=$" + par.costPerInf*dollarsPerInf)
    % 
    % 
    % % Decentralised problem
    % nexttile;
    % plot(t, resultsDecent(iScenario).R/sum(par.N))
    % hold on
    % set(gca, 'ColorOrderIndex', 1);
    % plot(t, results_u(iScenario).R/sum(par.N), '--')
    % plot(t, resultsDecent(iScenario).a, 'LineWidth', 2)
    % xlim([0 tHoriz])
    % grid on
    % if par.nGroups > 1
    %     lbls = ["R_"+groupLbls; "R_"+groupLbls+"_u"; "a_"+groupLbls ];
    % else
    %     lbls = ["R(t)", "R(t) (unmitigated)", "a(t)"];
    % end
    % legend(lbls, 'Location', 'southeast')
    % xlabel('time (days)')
    % title('(c) epidemic dynamics - decentralized control')
    % 
    % 
    % nexttile;
    % plot(t, (resultsDecent(iScenario).costInf + resultsDecent(iScenario).costCont)*dollarsPerInf/1e9)
    % hold on
    % plot(t, resultsDecent(iScenario).costInf*dollarsPerInf/1e9)
    % plot(t, resultsDecent(iScenario).costCont*dollarsPerInf/1e9)
    % set(gca, 'ColorOrderIndex', 1);
    % plot(t, (results_u(iScenario).costInf + results_u(iScenario).costCont)*dollarsPerInf/1e9, '--')
    % xlim([0 tHoriz])
    % grid on
    % if par.nGroups > 1
    %     lbls = ["cost_"+groupIDs; "cost_"+groupIDs+"_u"];
    % else
    %     lbls = ["total cost", "infection cost", "control cost", "total cost (unmitigated)"];
    % end
    % legend(lbls, 'Location', 'southeast')
    % ylabel('cumulative cost ($bn)')
    % xlabel('time (days)')
    % title('(d) costs - decentralized control')
    
    