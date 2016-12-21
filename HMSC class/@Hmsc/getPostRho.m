function res = getPostRho(m)

if ~isempty(m.truePar)
	valuesT = m.rhopw(m.truePar.rho, 1);
else
	valuesT = [];
end

values = nan(m.postSamN, 1 );
for i = 1:m.postSamN
	values(i,:) = m.rhopw(m.postSamVec(i).rho, 1);
end

res = {values, valuesT};

end