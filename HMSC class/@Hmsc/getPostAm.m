	function res = getPostAm(m)

if ~isempty(m.truePar)
	etasT = m.truePar.etas;
	lambdasT = m.truePar.lambdas;
	AmT = etasT*lambdasT;
	valuesT = AmT(:)';
else
	valuesT = [];
end

values = nan(m.postSamN, m.ncs*m.ns );
for j = 1:m.postSamN
	etas = m.postSamVec(j).etas;
	lambdas = m.postSamVec(j).lambdas;
	Am = etas*lambdas;
	values(j,:) = Am(:)';
end

res = {values, valuesT};

end