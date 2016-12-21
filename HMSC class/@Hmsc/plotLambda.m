function plotLambda(m, dispTrue, showToScreen, saveToFile, type)

for i=1:m.nr
	if m.factorCov(i)
		for k=1:m.ncr(i)
			res = m.getPostLambda(i, k);
			values = res{1};
			valuesT = res{2};
			label = strcat( 'lambda_', int2str(i), '_', int2str(k) );
			if strcmp(type, 'mix')
				m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
			end
			if strcmp(type, 'box')
				m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
			end
		end
	else
		res = m.getPostLambda(i, 0);
		values = res{1};
		valuesT = res{2};
		label = strcat( 'lambda_', int2str(i) );
		if strcmp(type, 'mix')
			m.plotMixing(values, valuesT, dispTrue, label, showToScreen, saveToFile)
		end
		if strcmp(type, 'box')
			m.plotBox(values, valuesT, dispTrue, label, showToScreen, saveToFile)
		end
	end
end


end