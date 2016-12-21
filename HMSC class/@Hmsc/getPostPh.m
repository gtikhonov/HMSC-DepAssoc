function res = getPostPh(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.ph';
else
	valuesT = [];
end

values = nan(m.postSamN, m.ns );
for j = 1:m.postSamN
	values(j,:) = m.postSamVec(j).ph';
end

res = {values, valuesT};

end