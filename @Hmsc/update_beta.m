function beta = update_beta(X,Xr,Xs,z,pi,ncr,eta,lambda,etas,lambdas,gamma,sigma,T,ph,iV,rho,phylogeny,iQg,detQg,outlierspecies,spatial,factorcov,speciesX,includeXs)
ns=size(T,1);
if speciesX
	nc=size(X{1},2);
else
	nc=size(X,2);
end
if length(spatial)==0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
	np = max(pi);
end
S = z;
if includeXs == 1
	S = S - Xs*etas*lambdas;
elseif includeXs == 2
	for j=1:ns
		S(:,j) = S(:,j) - Xs{j}*etas*lambdas(:,j);
	end
end

for r=1:nr
	eta1=eta{r};
	lambda1=lambda{r};
	if factorcov(r) > 0
		Xr1 = Xr{r};
		for k=1:ncr(r)
			if factorcov(r) == 1
				XrEta = diag(Xr1(:,k))*eta1;
				S = S-XrEta(pi(:,r),:)*lambda1(:,:,k);
			else
				for j=1:ns
					XrEta = diag(Xr1{j}(:,k))*eta1;
					S(:,j) = S(:,j) - XrEta(pi(:,r),:)*lambda1(:,j,k);
				end
			end
		end
	else
		S = S-eta1(pi(:,r),:)*lambda1;
	end
end
if ~phylogeny
	gaT = gamma'*T';
	for j=1:ns
		if speciesX
			X1 = X{j};
		else
			X1 = X;
		end
		Ube = inv(ph(j)*iV + (1/sigma(j,j))*(X1'*X1));
		mbe = Ube*(ph(j)*iV*gaT(:,j)+(1/sigma(j,j))*X1'*S(:,j));
		beta(:,j) = mvnrnd (mbe,Ube)';
	end
else
	if ~speciesX
		if outlierspecies
			PHI=(1./sqrt(ph))*(1./sqrt(ph))';
			Ube=inv(kron(X'*X,inv(sigma))+kron(eye(nc),1./PHI).*kron(iV,iQg(:,:,rho)));
			tmp1=(iQg(:,:,rho).*(1./PHI))*T*gamma*iV;
		else
			Ube=inv(kron(X'*X,inv(sigma))+kron(iV,iQg(:,:,rho)));
			tmp1=iQg(:,:,rho)*T*gamma*iV;
		end
		Ube=(Ube+Ube')/2;
		tmp2=inv(sigma)*S'*X;
		mbe=Ube*(tmp1(:)+tmp2(:));
		beta=mvnrnd(mbe,Ube)';
		beta=reshape(beta,ns,nc)';
	else
		tmp=X{1}'*X{1}/sigma(1,1);
		for j = 2:ns
			tmp = blkdiag(tmp,X{j}'*X{j}/sigma(j,j));
		end
		Ube=inv(tmp+kron(iQg(:,:,rho),iV));
		Ube=(Ube+Ube')/2;
		tmp1=iV*gamma'*T'*iQg(:,:,rho);
		tmp2 = zeros(nc,ns);
		for j=1:ns
			tmp2(:,j) = X{j}'*S(:,j)/sigma(j,j);
		end
		mbe=Ube*(tmp1(:)+tmp2(:));
		beta=mvnrnd(mbe,Ube)';
		beta=reshape(beta,nc,ns);
	end
end

