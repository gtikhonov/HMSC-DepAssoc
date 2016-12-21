function res = getPostAlpha(m, i)

if ~isempty(m.truePar)
	valuesT = m.alphapw{i}(m.truePar.alpha{i})';
else
	valuesT = [];
end

values = nan(m.postSamN, m.postSamVec(1).nf(i) );
for j = 1:m.postSamN
	values(j,:) = m.alphapw{i}(m.postSamVec(j).alpha{i})';
end

res = {values, valuesT};

end