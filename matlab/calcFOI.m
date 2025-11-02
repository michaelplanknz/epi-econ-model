function FOI = calcFOI(I, aFocal, aBG, par)

%prev = sum(I);
%TTI_effect = par.alpha_TTI + (1-par.alpha_TTI) * 1/(1 + exp(-(prev-par.TTI_max)/par.TTI_breadth));
TTI_effect = 1;

FOI = TTI_effect * ((aFocal .* par.Beta .* aBG') * I) ./ par.N;

 