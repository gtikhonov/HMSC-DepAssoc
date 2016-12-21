function res = getPostV(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.V(:);
else
	valuesT = [];
end

values = nan(m.postSamN, m.nc*m.nc );
for i = 1:m.postSamN
	values(i,:) = m.postSamVec(i).V(:);
end

res = {values, valuesT};

end