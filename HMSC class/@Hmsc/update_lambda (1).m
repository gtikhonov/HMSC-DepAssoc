function lambda =  update_lambda(X,Xr,Xs,z,beta,sigma,eta,lambda,etas,lambdas,delta,psijh,pi,nf,ncr,spatial,factorcov,speciesX,includeXs)
[ny ns] = size(z);
if length(spatial)==0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
	np = max(pi);
end
for r1=1:nr
	S = z;
	if speciesX
		for j=1:ns
			S(:,j) = S(:,j) - X{j}*beta(:,j);
		end
	else
		S = S - X*beta;
	end
	if includeXs == 1
		S = S - Xs*etas*lambdas;
	elseif includeXs == 2
		for j=1:ns
			S(:,j) = S(:,j) - Xs{j}*etas*lambdas(:,j);
		end
	end
	for r2=1:nr
		if ~(r1==r2)  
			eta1=eta{r2};
			lambda1=lambda{r2};
			if factorcov(r2) > 0
				Xr1 = Xr{r2};
				for k=1:ncr(r2)
					if factorcov(r2) == 1
						XrEta = repmat(Xr1(:,k), 1, nf(r2)).*eta1;
						S = S-XrEta(pi(:,r2),:)*lambda1(:,:,k);
					else
						for j=1:ns
							XrEta = repmat(Xr1{j}(:,k), 1, nf(r2)).*eta1;
							S(:,j) = S(:,j) - XrEta(pi(:,r2),:)*lambda1(:,j,k);
						end
					end
				end
			else
				S = S-eta1(pi(:,r2),:)*lambda1;
			end
		end
	end
	lpi = pi(:,r1);
	eta1=eta{r1};
	lambda1=lambda{r1};
	delta1=delta{r1};
	psijh1=psijh{r1};
	tauh = cumprod(delta1);
	if factorcov(r1) == 1
		Xr1 = Xr{r1};
		Plam = zeros( ns, nf(r1), ncr(r1) );
		XrEta=[];
		for k=1:ncr(r1)
			XrEta = [ XrEta, repmat(Xr1(:,k), 1, nf(r1)).*eta1 ];
			Plam(:,:,k) = bsxfun(@times,psijh1(:,:,k),tauh(:,k)');
		end
		eta1 = XrEta;
		eta1=eta1(lpi,:);
		eta2 = eta1'*eta1;
		for j = 1:ns
			Qlam = diag(Plam(j,:)) + (1/sigma(j,j))*eta2;
			blam = (1/sigma(j,j))*eta1'*S(:,j);
			Llam = chol(Qlam,'lower');
			zlam = normrnd(0,1,nf(r1)*ncr(r1),1);
			vlam = Llam\blam;
			mlam = Llam'\vlam;
			ylam = Llam'\zlam;
			lambda1(:,j,:) = reshape((ylam + mlam), [nf(r1),1,ncr(r1)]);
		end
	elseif factorcov(r1) == 2
		eta0 = eta1;
		Xr1 = Xr{r1};
		Plam = zeros( ns, nf(r1), ncr(r1) );
		for k=1:ncr(r1)
			Plam(:,:,k) = bsxfun(@times,psijh1(:,:,k),tauh(:,k)');
		end
		for j = 1:ns
			XrEta=[];
			for k=1:ncr(r1)
				XrEta = [ XrEta, repmat(Xr1{j}(:,k), 1, nf(r1)).*eta0 ];
			end
			eta1 = XrEta;
			eta1=eta1(lpi,:);
			eta2 = eta1'*eta1;
			Qlam = diag(Plam(j,:)) + (1/sigma(j,j))*eta2;
			blam = (1/sigma(j,j))*eta1'*S(:,j);
			Llam = chol(Qlam,'lower');
			zlam = normrnd(0,1,nf(r1)*ncr(r1),1);
			vlam = Llam\blam;
			mlam = Llam'\vlam;
			ylam = Llam'\zlam;
			lambda1(:,j,:) = reshape((ylam + mlam), [nf(r1),1,ncr(r1)]);
		end
	else
		Plam = bsxfun(@times,psijh1,tauh');
		eta1=eta1(lpi,:);
		eta2 = eta1'*eta1;
		for j = 1:ns
			Qlam = diag(Plam(j,:)) + (1/sigma(j,j))*eta2;
			blam = (1/sigma(j,j))*eta1'*S(:,j);
			Llam = chol(Qlam,'lower');
			zlam = normrnd(0,1,nf(r1),1);
			vlam = Llam\blam; mlam = Llam'\vlam; ylam = Llam'\zlam;
			lambda1(:,j) = (ylam + mlam);
		end
	end
	lambda{r1}=lambda1;
end