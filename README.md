# Joint economic and epidemiological modelling of alternative pandemic response strategies

This repository contains Matlab code to reproduce the results in the article [Joint economic and epidemiological modelling of alternative pandemic response strategies](https://arxiv.org/abs/2512.08355).

# Abstract

In an emerging pandemic, policymakers need to make decisions with limited information and requiring trade-offs between the health impact of the pandemic and the economic costs of the response. Most mathematical models have focused on direct health impacts, neglecting the economic costs of control measures. Here, we introduce a framework that captures both health and economic costs and compare the expected aggregate costs of alternative strategies across a range of epidemiological and economic parameters. We find that for diseases with low severity, mitigation tends to be the most cost-effective option. For more severe diseases, suppression tends to be most cost effective if the basic reproduction number R_0 is relatively low, while elimination tends to be more cost-effective if R_0 is high. We use the example of New Zealand's response to the Covid-19 pandemic in 2020 to anchor our framework to a real-world case study. We find that parameter estimates for Covid-19 in New Zealand put it close to or above the threshold at which elimination becomes more cost-effective than mitigation. We conclude that our proposed framework holds promise as a decision-support tool for pandemic threats, with further work needed to account for population heterogeneity and other factors relevant to decision-making. 




# Version history

* 9 December 2025 [pre-print](https://arxiv.org/abs/2512.08355v1) - results generated using the version of this repository tagged `v1.0`. 

# How to use this repository

* Navigate to the `Matlab/` folder and run the top-level script `main.m`. This will run the model across a range of parameter combinations and save the results in the `/output/` folder as `output/results.mat`.
* Run the script `postProcess.m`. This will read in the previously saved results and calculate the suppression and elimination costs for the corresponding parameter range, and create various arrays for plotting. Results will be saved as `output/results_for_plots.mat`.
* Run the script `plotGraphs.m`. This will read the previously saved results back in and plot the graphs. If the variable `saveFlag` is set to true, the Figures will be saved as .png files in the `figures/` folder.
* If you want to change any parameter values, edit the function `getPar()` which sets all model parameters. 


