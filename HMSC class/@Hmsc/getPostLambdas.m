function res = getPostLambdas(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.lambdas(:)';
else
	valuesT = [];
end

values = nan(m.postSamN, m.postSamVec(1).nfs*m.ns );
for j = 1:m.postSamN
	values(j,:) = m.postSamVec(j).lambdas(:)';
end

res = {values, valuesT};

end