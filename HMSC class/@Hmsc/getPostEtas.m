function res = getPostEtas(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.etas(:)';
else
	valuesT = [];
end

values = nan(m.postSamN, m.ncs*m.postSamVec(1).nfs );
for j = 1:m.postSamN
	values(j,:) = m.postSamVec(j).etas(:)';
end

res = {values, valuesT};

end