function res = getPostEta(m, i)

if ~isempty(m.truePar)
	valuesT = m.truePar.eta{i}(:)';
else
	valuesT = [];
end

% not implemented - different number of factors....

values = nan(m.postSamN, m.np(i) * m.postSamVec(1).nf(i) );
for j = 1:m.postSamN
	values(j,:) = m.postSamVec(j).eta{i}(:)';
end

res = {values, valuesT};

end