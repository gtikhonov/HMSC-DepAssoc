function genY(m, misFraction)
p = m.truePar;

if m.speciesX
	Ez = zeros(m.ny, m.ns);
	for r = 1:m.ns
		Ez(:,r) = m.X{r}*p.beta(:,r);
	end
else
	Ez = m.X*p.beta;
end

for r = 1:m.nr
	eta1 = p.eta{r};
	lambda1 = p.lambda{r};
	if m.factorCov(r) == 1
		for k = 1:m.ncr(r)
			Xreta = repmat(m.Xr{r}(:,k), 1, p.nf(r)).*eta1;
			Ez = Ez + Xreta(m.pi(:,r),:)*lambda1(:,:,k);
		end
	elseif m.factorCov(r) == 2
		for j = 1:m.ns
			for k = 1:m.ncr(r)
				Xreta = repmat(m.Xr{r}{j}(:,k), 1, p.nf(r)).*eta1;
				Ez(:,j) = Ez(:,j) + Xreta(m.pi(:,r),:)*lambda1(:,j,k);
			end
		end
	else
		Ez = Ez + eta1(m.pi(:,r),:)*lambda1;
	end
end

if m.includeXs == 1
	Xel = m.Xs*p.etas*p.lambdas;
	Ez = Ez+Xel;
elseif m.includeXs == 2
	for r=1:m.ns
		Ez(:,r) = Ez(:,r) + m.Xs{r}*p.etas*p.lambdas(:,r);
	end
end

eps = zeros(m.ny, m.ns);
for r = 1:m.ny
	eps(r,:) = mvnrnd(zeros(1,m.ns), p.sigma);
end
k=ones(m.ny, m.ns);
for r = 1:m.ns
	if m.dist(r,3) == 1
		k(:,r) = max(Ez(:,r),1).^m.dist(r,4);
	end
	if m.dist(r,3) == 2
		k(:,r) = exp(Ez(:,r)).^m.dist(r,4);
	end
end
z = Ez + k.*eps;
Y = z;
for j = 1:m.ns
	if(m.dist(j,1) == 2)
		Y(:,j) = z(:,j)>0;
	end
	if(m.dist(j,1) == 3)
		Y(:,j) = poissrnd(exp(z(:,j)));
	end
	if(m.dist(j,1) == 4)
		Y(:,j) = poissrnd(max(0,z(:,j)));
	end
	for r = 1:m.ny
		if rand < misFraction
			Y(r,j) = NaN;
		end
	end
end
m.setY(Y);
end