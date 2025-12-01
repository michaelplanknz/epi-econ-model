clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in prevously saved results from optimization and calculates cost of
% mitigation and suppression and calculates various quantities for plotting

% Set output folder location
outFolder = '../output/';


% Load previously saved model results
fIn = outFolder + "results.mat";
load(fIn);

% Values for final heat map plot
Beta_fix = 0.6;
costPerInf_fix = 0.6;
tDetRet_vals = 0:0.005:0.2;
alpha_TTI_vals = 0.5:0.01:1;


% Upper limit for y axis for cost plots for each beta
cUpper = [40 80];

% Main results use cost per infection measured in $10,000s
% This scaling factor converts all costs to $
dollarsPerInf = 10000;


% Time horizon to use for coast comparisons and plots (need to be <= par.tMax
tHoriz = 600;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate results for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialise matrices for heat map results
costUnmit = zeros(size(Beta_mat));
costDecent = zeros(size(Beta_mat));
costMit = zeros(size(Beta_mat));
costSup = zeros(size(Beta_mat));
costElim = zeros(size(Beta_mat));
stratCode = nan(size(Beta_mat));
tCrit = nan(size(Beta_mat));
costInfElim = nan(size(Beta_mat));

% Calculate strategy costs across scenarios
parInd = par;
for iScenario = 1:nScenarios

    %Set scneario-dependent parameters
    parInd.Beta = Beta_list(iScenario);
    parInd.costPerInf = costPerInf_list(iScenario);

    % Get matrix indices for this iScenario
    [iRow, jCol] = ind2sub(size(Beta_mat), iScenario);

    % Fill matrix entries
    costUnmit(iRow, jCol) = results_u(iScenario).costInf(end);
    costDecent(iRow, jCol) = resultsDecent(iScenario).costInf(end) + resultsDecent(iScenario).costCont(end);
    costMit(iRow, jCol) = resultsCent(iScenario).costInf(end) + resultsCent(iScenario).costCont(end);

     % Compute elimination costs
    [CElimRate, aOptElim, ~, ~, CSupRate, ~, incRateElim] = calcElimCost(t, parInd);
    costElim(iRow, jCol) = CElimRate*tHoriz;
    costSup(iRow, jCol) = CSupRate*tHoriz;
    
    % Calculate the cost associated with infections from border-related
    % outbreaks under the elimination strategy (this is ignored in
    % subsequent calculations, but this is to check that it is negligible)
    % [NB only calculate this if R0*alpha > 1, otherwise no reduction in
    % activity is needed and suppression will always be preferred over
    % elimination]
    if parInd.Beta/parInd.Gamma*parInd.alpha_TTI > 1 
        % Calculation is incidence rate (infections per unit time) * time
        % horizon * cost per infection
        costInfElim(iRow, jCol) = parInd.costPerInf * incRateElim*tHoriz;
    end

    % Record code for optimal strategy: 1 = mitigation, 2 = suppression, 3 = elimination
    [~, stratCode(iRow, jCol)] = min([costMit(iRow, jCol), costSup(iRow, jCol), costElim(iRow, jCol)]);

    % Set code = 0 when mitigation response does not appreciably reduce
    % activity or cost is less than 5% below decentralised cost
    if stratCode(iRow, jCol) == 1 & min(resultsCent(iScenario).a) >= 0.999 %| costMit(iRow, jCol) >= 0.95*costDecent(iRow, jCol))
        stratCode(iRow, jCol) = 0;
    end




    % Calculate threshold time for elimination
    tCrit(iRow, jCol) = costMit(iRow, jCol)/min(CElimRate, CSupRate);

end


% For fixed (R0, k) calculate tCrit across detect-return ratio and TII effectiveness
parInd.Beta = Beta_fix;
parInd.costPerInf = costPerInf_fix;

% Get mitigation cost for comparison
iBeta = find(abs(Beta_vals - Beta_fix) < 1e-8);
iCost = find(abs(costPerInf_vals - costPerInf_fix) < 1e-8);
costComp = costMit(iCost, iBeta);

% Loop over tDetRet and alpha_TTI
nRows = length(alpha_TTI_vals);
mCols = length(tDetRet_vals);
tCrit2 = zeros(nRows, mCols);
for iRow = 1:nRows
    parInd.alpha_TTI = alpha_TTI_vals(iRow);
    for jCol = 1:mCols
        parInd.r = tDetRet_vals(jCol)/parInd.tDet;
        [CElimRate, ~, ~, ~, CSupRate] = calcElimCost(t, parInd);    
        tCrit2(iRow, jCol) = costComp/min(CElimRate, CSupRate);
    end
end



% Save processed results
fOut = outFolder + "results_for_plots.mat";
save(fOut);

