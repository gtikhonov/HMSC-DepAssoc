function plotGamma(m, dispTrue, showToScreen, saveToFile, type)

res = m.getPostGamma();
values = res{1};
valuesT = res{2};

label = 'gamma';
if strcmp(type, 'mix')
	m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end
if strcmp(type, 'box')
	m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
end

end