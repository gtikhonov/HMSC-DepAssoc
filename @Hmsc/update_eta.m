function eta = update_eta(z,X,Xr,Xs,beta,sigma,eta,alpha,lambda,etas,lambdas,nf,pi,ncr,spatial,factorcov,iWg,speciesX,includeXs)
[ny ns] = size(z);
isigma = inv(sigma);
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
						XrEta = diag(Xr1(:,k))*eta1;
						S = S-XrEta(pi(:,r2),:)*lambda1(:,:,k);
					else
						for j=1:ns
							XrEta = diag(Xr1{j}(:,k))*eta1;
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
	lambda1=lambda{r1};
	if factorcov(r1)
		Xr1 = Xr{r1};
		if ~spatial(r1)
			eta1=eta{r1};
			if sum(lpi==(1:ny)')==ny % i1 corresponds to residuals, use a faster updater
				for p = 1:ny % change to parfor for time demanding works
					lambdaX = zeros(nf(r1), ns);
					for k = 1:ncr(r1)
						if factorcov(r1) == 1
							lambdaX = lambdaX + Xr1(p,k)*lambda1(:,:,k);
						else
							for j=1:ns
								lambdaX(:,j) = lambdaX(:,j) + Xr1{j}(p,k)*lambda1(:,j,k);
							end
						end
					end
					lambda2 = (lambdaX*isigma*lambdaX');
					Veta1 = eye(nf(r1)) + lambda2;
					Tx = cholcov(Veta1); [Q,R] = qr(Tx);
					iR = inv(R); Veta = iR*iR';  % Veta = inv(Veta1)
					Meta = S(p,:)*isigma*lambdaX'*Veta;
					eta1(p,:) = Meta + normrnd(0,1,[1,nf(r1)])*iR';
				end
			else
				for p = 1:np(r1)
					lambdaX = zeros(nf(r1), ns);
					for k = 1:ncr(r1)
						if factorcov(r1) == 1
							lambdaX = lambdaX + Xr1(p,k)*lambda1(:,:,k);
						else
							for j=1:ns
								lambdaX(:,j) = lambdaX(:,j) + Xr1{j}(p,k)*lambda1(:,j,k);
							end
						end
					end
					lambda2 = (lambdaX*isigma*lambdaX');
					rows = lpi==p;
					Veta1 = eye(nf(r1)) + lambda2*sum(rows);
					Tx = cholcov(Veta1);
					[Q,R] = qr(Tx);
					iR = inv(R); Veta = iR*iR';  % Veta = inv(Veta1)
					Meta = sum(S(rows,:),1)*isigma*lambdaX'*Veta;
					eta1(p,:) = Meta + normrnd(0,1,[1,nf(r1)])*iR';
				end
			end			
		else % SPATIAL LATENT FACTORS
			iWg1 = iWg{r1};
			alpha1=alpha{r1};
			iWs = zeros( np(r1)*nf(r1) );
			for h=1:nf(r1)
				iWs((h-1)*np(r1)+(1:np(r1)), (h-1)*np(r1)+(1:np(r1))) = iWg1(:,:,alpha1(h));
			end
			Xrk = cell(1, ncr(r1));
			isigmaLambdak = cell(1, ncr(r1));
			for k = 1:ncr(r1)
				if factorcov(r1) == 1
					Xrk1 = diag(Xr1(:,k));
					Xrk{k} = Xrk1(lpi,:);
				else
					fprintf('Not implemented yet\n');
					Xrk1 = cell(1, ns);
					for j = 1:ns
						Xrk1{j} = diag(Xr1{j}(:,k));
						Xrk1{j} = Xrk1{j}(lpi,:);
					end
					Xrk{k} = Xrk1;
				end
				isigmaLambdak{k} = isigma*lambda1(:,:,k)';
			end
			Ueta1 = iWs;
			for k1 = 1:ncr(r1)
				Xrk1 = Xrk{k1};
				for k2 = 1:ncr(r1)
					Xrk2 = Xrk{k2};
					if factorcov(r1) == 1
						Ueta1 = Ueta1 + kron(lambda1(:,:,k1)*isigmaLambdak{k2}, Xrk1'*Xrk2);
					else
						fprintf('Not implemented yet\n');
					end
				end
			end
			Ueta = inv(Ueta1);
			Ueta = (Ueta+Ueta')/2;
			Meta = zeros(nf(r1)*np(r1), 1);
			for k1 = 1:ncr(r1)
				Xrk1 = Xrk{k1};
				vec = Xrk1' * S * isigmaLambdak{k1};
				Meta = Meta + Ueta*vec(:);
			end;
			feta = mvnrnd(Meta,Ueta);
			eta1 = reshape(feta,[np(r1),nf(r1)]);
		end
	else
		if spatial(r1)==0 % NON SPATIAL LATENT FACTORS
			if length(unique(lpi))==ny % i1 corresponds to residuals, use a faster updater 
				% if sum(lpi==(1:ny)')==ny
				% 2016.07.30 - this change is rather experimental, as I cannot clearly predic whether it can lead to
				% misbehaviour in some situation or not. However it should work correctly for estimation and speeds up the
				% coditional predictions a lot.
				Veta1 = eye(nf(r1)) + lambda1*isigma*lambda1';
				Tx = cholcov(Veta1);
				[Q,R] = qr(Tx);
				iR = inv(R);
				Veta = iR*iR';
				Meta = S*isigma;
				Meta = Meta*lambda1';
				Meta = Meta*Veta;                        % ny x k
				eta1(lpi,:) = Meta + normrnd(0,1,[ny,nf(r1)])*iR';       % update eta1 in a block
			else % i1 corresponds to a random effect, use a more general updater
				eta1=eta{r1};
				lambda2 = (lambda1*isigma*lambda1');
				for p = 1:np(r1)
					rows = lpi==p;
					Veta1 = eye(nf(r1)) + lambda2*sum(rows);
					Tx = cholcov(Veta1); 	[Q,R] = qr(Tx);
					iR = inv(R); Veta = iR*iR';   % Veta = inv(Veta1)
					Meta = sum(S(rows,:),1)*isigma*lambda1'*Veta;
					eta1(p,:) = Meta + normrnd(0,1,[1,nf(r1)])*iR';
				end
			end
		else % SPATIAL LATENT FACTORS
			iWg1 = iWg{r1};
			alpha1=alpha{r1};
			iWs = zeros( np(r1)*nf(r1) );
			for h=1:nf(r1)
				iWs((h-1)*np(r1)+(1:np(r1)), (h-1)*np(r1)+(1:np(r1))) = iWg1(:,:,alpha1(h));
			end
			if sum(lpi==(1:ny)')==ny % i1 corresponds to spatial residuals, use a faster updater
				tmp = isigma*lambda1';
				tmp1 = lambda1*tmp;
				tmp1s = kron(tmp1,eye(ny));
				iUeta = iWs + tmp1s;
				Ueta = inv(iUeta);
				Ueta = (Ueta+Ueta')/2;
				fS = S*tmp;
				Meta = Ueta*fS(:);
				feta = mvnrnd(Meta,Ueta);
				eta1 = reshape(feta,[np(r1),nf(r1)]);
			else % i1 corresponds to a spatial random effect, use a more general updater
				Pi = eye(np(r1));
				Pi = Pi(lpi, :);
				tmp = isigma*lambda1';
				tmp1 = lambda1*tmp;
				tmp1s = kron(tmp1,Pi'*Pi);
				iUeta = iWs + tmp1s;
				Ueta = inv(iUeta);
				Ueta = (Ueta+Ueta')/2;
				fS = Pi'*S*tmp;
				Meta = Ueta*fS(:);
				feta = mvnrnd(Meta,Ueta);
				eta1 = reshape(feta,[np(r1),nf(r1)]);
			end
		end
	end
	eta{r1} = eta1;
end
