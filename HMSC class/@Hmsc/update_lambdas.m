function lambdas=update_lambdas(X,Xs,Xr,pi,ncr,z,lambdas,etas,beta,sigma,eta,lambda,psijhs,deltas,spatial,speciesX,includeXs,factorcov)
[ny ncs]=size(Xs);
[nfs ns]=size(lambdas);
if length(spatial)==0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
	np = max(pi);
end
S = z;
if speciesX
	for i=1:ns
		S(:,i) = S(:,i) - X{i}*beta(:,i);
	end
else
	S = S - X*beta;
end
for i1=1:nr
	eta1=eta{i1};
	lambda1=lambda{i1};
	if factorcov(i1)
		Xr1 = Xr{i1};
		for k1=1:ncr(i1)
			XrEta = repmat(Xr1(:,k), 1, nf(r)).*eta1;
			S = S-XrEta(pi(:,i1),:)*lambda1(:,:,k1);
		end
	else
		S = S-eta1(pi(:,i1),:)*lambda1;
	end
end

taus = cumprod(deltas);
if includeXs == 2
	for j=1:ns
		Dj=diag(1./(psijhs(j,:).*taus'));
		tmp=Xs{j}*etas;
		Uj=inv(inv(Dj)+(1/sigma(j,j))*(tmp'*tmp));
		mj=Uj*((1/sigma(j,j))*tmp'*S(:,j));
		lambdas(:,j)=mvnrnd(mj,Uj);
	end
else
	tmp=Xs*etas;
	for j=1:ns
		Dj=diag(1./(psijhs(j,:).*taus'));
		Uj=inv(inv(Dj)+(1/sigma(j,j))*(tmp'*tmp));
		mj=Uj*((1/sigma(j,j))*tmp'*S(:,j));
		lambdas(:,j)=mvnrnd(mj,Uj);
	end
end
