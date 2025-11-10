function FOI = calcFOI(I, aFocal, aBG, par)


FOI = ((aFocal .* par.Beta .* aBG') * I) ./ par.N;

 