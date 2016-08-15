function genY(m, misFraction)
p = m.truePar;

if m.speciesX
	Ez = zeros(m.ny, m.ns);
	for i = 1:m.ns
		Ez(:,i) = m.X{i}*p.beta(:,i);
	end
else
	Ez = m.X*p.beta;
end

for i = 1:m.nr
	eta1 = p.eta{i};
	lambda1 = p.lambda{i};
	if m.factorCov(i) == 1
		for k = 1:m.ncr(i)
			Xreta = diag(m.Xr{i}(:,k))*eta1;
			Ez = Ez + Xreta(m.pi(:,i),:)*lambda1(:,:,k);
		end
	elseif m.factorCov(i) == 2
		for j = 1:m.ns
			for k = 1:m.ncr(i)
				Xreta = diag(m.Xr{i}{j}(:,k))*eta1;
				Ez(:,j) = Ez(:,j) + Xreta(m.pi(:,i),:)*lambda1(:,j,k);
			end
		end
	else
		Ez = Ez + eta1(m.pi(:,i),:)*lambda1;
	end
end

if m.includeXs == 1
	Xel = m.Xs*p.etas*p.lambdas;
	Ez = Ez+Xel;
elseif m.includeXs == 2
	for i=1:m.ns
		Ez(:,i) = Ez(:,i) + m.Xs{i}*p.etas*p.lambdas(:,i);
	end
end

eps = zeros(m.ny, m.ns);
for i = 1:m.ny
	eps(i,:) = mvnrnd(zeros(1,m.ns), p.sigma);
end
k=ones(m.ny, m.ns);
for i = 1:m.ns
	if m.dist(i,3) == 1
		k(:,i) = max(Ez(:,i),1).^m.dist(i,4);
	end
	if m.dist(i,3) == 2
		k(:,i) = exp(Ez(:,i)).^m.dist(i,4);
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
	for i = 1:m.ny
		if rand < misFraction
			Y(i,j) = NaN;
		end
	end
end
m.setY(Y);
end