function sigma = update_sigma(X,Xr,Xs,sigma,beta,eta,lambda,etas,lambdas,nf,pi,ncr,dist,asigma,bsigma,z,spatial,factorcov,speciesX,includeXs)
if length(spatial)==0
	nr = 0;
else
	[ny nr]=size(pi);
end
[ny ns] = size(z);
if speciesX
	Ez = zeros(ny,ns);
	for j=1:ns
		Ez(:,j)=X{j}*beta(:,j);
	end
else
	Ez = X*beta;
end

if includeXs == 1
	Ez = Ez + Xs*etas*lambdas;
elseif includeXs == 2
	for i=j:ns
		Ez(:,j) = Ez(:,j) + Xs{j}*etas*lambdas(:,j);
	end
end

for r=1:nr
	eta1=eta{r};
	lambda1=lambda{r};
	if factorcov(r) > 0
		Xr1 = Xr{r};
		for k=1:ncr(r)
			if factorcov(r) == 1
				XrEta = repmat(Xr1(:,k), 1, nf(r)).*eta1;
				Ez = Ez+XrEta(pi(:,r),:)*lambda1(:,:,k);
			else
				for j=1:ns
					XrEta = repmat(Xr1{j}(:,k), 1, nf(r)).*eta1;
					Ez(:,j) = Ez(:,j) + XrEta(pi(:,r),:)*lambda1(:,j,k);
				end
			end
		end
	else
		Ez = Ez+eta1(pi(:,r),:)*lambda1;
	end
end
eps=z-Ez;
for j=1:ns
	if (dist(j,2))
		tmp=1/gamrnd(asigma(j)+ny/2,1./(bsigma(j)+sum(eps(:,j).^2)/2));
		if ~isnan(tmp)
			sigma(j,j)=tmp;
		end
	end
end
