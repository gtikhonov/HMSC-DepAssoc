function res = getPostGamma(m)

if ~isempty(m.truePar)
	valuesT = m.truePar.gamma(:)';
else
	valuesT = [];
end

values = nan(m.postSamN, m.nt*m.nc );
for i = 1:m.postSamN
	values(i,:) = m.postSamVec(i).gamma(:)';
end

res = {values, valuesT};

end