function plotAlpha(m, dispTrue, showToScreen, saveToFile, type)

for i=1:m.nr
	if m.spatial(i)
		res = m.getPostAlpha(i);
		values = res{1};
		valuesT = res{2};
		label = strcat( 'alpha_', int2str(i) );
		if strcmp(type, 'mix')
			m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
		end
		if strcmp(type, 'box')
			m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
		end
	end
end


end