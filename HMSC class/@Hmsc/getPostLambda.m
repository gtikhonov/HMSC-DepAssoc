function res = getPostLambda(m, i, k)

if ~isempty(m.truePar)
	lambdaT = m.truePar.lambda;
else
	lambdaT = [];
end

values = nan(m.postSamN, m.postSamVec(1).nf(i)*m.ns );

if m.factorCov(i)
	valuesT=[];
	if ~isempty(m.truePar)
		valuesT=lambdaT{i}(:,:,k);
		valuesT = valuesT(:)';
	end
	for j = 1:m.postSamN
		lambda = m.postSamVec(j).lambda{i}(:,:,k);
		values(j,:) = lambda(:)';
	end
else
	valuesT=[];
	if ~isempty(m.truePar)
		valuesT=lambdaT{i}(:)';
	end
	for j = 1:m.postSamN
		lambda = m.postSamVec(j).lambda{i};
		values(j,:) = lambda(:)';
	end
end

res = {values, valuesT};

end