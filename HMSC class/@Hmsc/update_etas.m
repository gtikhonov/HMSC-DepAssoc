function etas= update_etas(X,Xs,Xr,pi,ncr,z,lambdas,beta,sigma,eta,lambda,spatial,speciesX,includeXs,factorcov)
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
for r1=1:nr
	eta1=eta{r1};
	lambda1=lambda{r1};
	if factorcov(r1)
		Xr1 = Xr{r1};
		for k1=1:ncr(r1)
			XrEta = repmat(Xr1(:,k), 1, nf(r1)).*eta1;
			S = S-XrEta(pi(:,r1),:)*lambda1(:,:,k1);
		end
	else
		S = S-eta1(pi(:,r1),:)*lambda1;
	end
end

isigma=inv(sigma);
if includeXs == 2
	[ny, ncs]=size(Xs{1});
	Ksi = nan(ny*ns, nfs*ncs);
	tmp = Ksi';
	for i = 1:ns
		rows = (i-1)*ny+(1:ny);
		Ksi(rows,:) = kron(lambdas(:,i)', Xs{i});
		tmp(:,rows) = isigma(i,i)*(Ksi(rows,:))';
	end
	% 	tmp = Ksi'*kron(isigma, eye(ny));
	Uetas = inv( eye(ncs*nfs)+tmp*Ksi );
	Uetas = (Uetas+Uetas')/2;
	metas = Uetas*tmp*S(:);
else
	[ny, ncs]=size(Xs);
	Uetas=inv( eye(ncs*nfs)+kron(lambdas*isigma*lambdas',Xs'*Xs) );
	Uetas=(Uetas+Uetas')/2;
	tmp=Xs'*S*isigma*lambdas';
	metas=Uetas*(tmp(:));
end
tmp = mvnrnd(metas, Uetas);
etas=reshape(tmp,[ncs,nfs]);