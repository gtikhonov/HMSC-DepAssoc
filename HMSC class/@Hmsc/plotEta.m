function plotEta(m, dispTrue, showToScreen, saveToFile, type)

for k = 1:m.nr
	res = m.getPostEta(k);
	values = res{1};
	valuesT = res{2};
	
	label = strcat('eta_', int2str(k));
	if strcmp(type, 'mix')
		m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
	end
	if strcmp(type, 'box')
		m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
	end
end;

end