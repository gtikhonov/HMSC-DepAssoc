function plotAll(m, dispTrue, type)

m.plotSigma(dispTrue, false, true, type)
m.plotGamma(dispTrue, false, true, type)
m.plotBeta(dispTrue, false, true, type)
m.plotEta(dispTrue, false, true, type)
m.plotLambda(dispTrue, false, true, type)
if m.phylogeny
	m.plotRho(dispTrue, false, true, type)
end
if any(m.spatial)
	m.plotAlpha(dispTrue, false, true, type)
end
if m.outlierSpecies
	m.plotPh(dispTrue, false, true, type)
end
if m.includeXs
	m.plotLambdas(dispTrue, false, true, type)
	m.plotEtas(dispTrue, false, true, type)
	m.plotAm(dispTrue, false, true, type)
end

end