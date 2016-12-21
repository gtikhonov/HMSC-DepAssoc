function res = getPostSigma(m)

if ~isempty(m.truePar)
	valuesT = diag(m.truePar.sigma);
else
	valuesT = [];
end

values = nan(m.postSamN, m.ns );
for i = 1:m.postSamN
	values(i,:) = diag(m.postSamVec(i).sigma);
end

res = {values, valuesT};

end