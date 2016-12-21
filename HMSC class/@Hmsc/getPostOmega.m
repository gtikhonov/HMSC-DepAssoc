function res = getPostOmega(m, i, x)

if ~isempty(m.truePar)
	lambdaT = m.truePar.lambda;
else
	lambdaT = [];
end

values = nan(m.postSamN, m.ns*m.ns);

if m.factorCov(i)
	if any(size(x) ~= [1, m.ncr(i)])
		error('HMSC: the length of x vector is not equal to number of covariates')
	end
	valuesT=[];
	if ~isempty(m.truePar)
		lam = zeros( size(lambdaT{i}(:,:,1)) );
		for k = 1:m.ncr(i)
			lam = lam + x(k)*lambdaT{i}(:,:,k);
		end
		RR = lam'*lam;
		valuesT = RR(:)';
	end
	for j = 1:m.postSamN
		lambda = m.postSamVec(j).lambda{i};
		lam = zeros( size(lambda(:,:,1)) );
		for k = 1:m.ncr(i)
			lam = lam + x(k)*lambda(:,:,k);
		end
		RR = lam'*lam;
		values(j,:) = RR(:)';
	end
else
	if ~isequal(x, [])
		error('HMSC: non-empty x vector was given to latent factor level without covariates')
	end
	valuesT=[];
	if ~isempty(m.truePar)
		RR=lambdaT{i}'*lambdaT{i};
		fR = RR(:)';
		valuesT=fR;
	end
	for j=1:m.postSamN
		lambda=m.postSamVec(j).lambda{i};
		RR=lambda'*lambda;
		values(j,:) = RR(:)';
	end
end

res = {values, valuesT};

end