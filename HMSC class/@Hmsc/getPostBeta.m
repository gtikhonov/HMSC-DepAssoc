function res = getPostBeta(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.beta(:)';
else
	valuesT = [];
end

values = nan(m.postSamN, m.nc*m.ns );
for i = 1:m.postSamN
	values(i,:) = m.postSamVec(i).beta(:)';
end

res = {values, valuesT};

end