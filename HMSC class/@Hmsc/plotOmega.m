function plotOmega(m, Xr, maxPairs, dispTrue, showToScreen, saveToFile, type)

for r=1:m.nr
	res = m.getPostOmega(r, Xr{r});
	values = res{1};
	valuesT = res{2};
	
	if size(values, 2) > maxPairs
		rng(0);
		ind = sort(randsample(size(values, 2), maxPairs));
		values = values(:,ind);
		if ~isempty(valuesT)
			valuesT = valuesT(:,ind);
		end
	end
	
	label = ['omega_', int2str(r)];
	if strcmp(type, 'mix')
		m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
	end
	if strcmp(type, 'box')
		m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
	end
end

end