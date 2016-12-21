function plotRho(m, dispTrue, showToScreen, saveToFile, type)

res = m.getPostRho();
values = res{1};
valuesT = res{2};

label = 'rho';
if strcmp(type, 'mix')
	m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end
if strcmp(type, 'box')
	m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end

end