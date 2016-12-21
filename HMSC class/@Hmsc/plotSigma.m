function plotSigma(m, dispTrue, showToScreen, saveToFile, type)

res = m.getPostSigma();
values = res{1};
valuesT = res{2};

label = 'sigma';
if strcmp(type, 'mix')
	m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end
if strcmp(type, 'box')
	m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end

end