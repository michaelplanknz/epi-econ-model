# Joint economic and epidemiological modelling of alternative pandemic response strategies

This repository contains Matlab code to reproduce the results in the article ``Joint economic and epidemiological modelling of alternative pandemic response strategies''.

# How to use this repository

* Navigate to the `Matlab/` folder and run the top-level script `main.m`. This will run the model across a range of parameter combinations and save the results in the `output/` folder .
* Run the script `plotGraphs.m`. This will read the previously saved results back in and plot the graphs. If the variable `saveFlag` is set to true, the Figures will be saved as .png files in the `figures/` folder.
* If you want to change any parameter values, edit the function `getPar()` which sets all model parameters. 


