function [X, covScale, Xs, sCovScale, Xr, factorCovScale, beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs] = computeInitialValues(m, start)

[covScale, X] = Hmsc.scaleMatrix(m.X, m.covScaleFlag, m.speciesX);
if m.includeXs
	[sCovScale, Xs] = Hmsc.scaleMatrix(m.Xs, m.sCovScaleFlag, m.includeXs-1);
else
	sCovScale = [];
	Xs = [];
end
factorCovScale = cell(1, m.nr);
Xr = cell(1, m.nr);
for r=1:m.nr
	if m.factorCov(r)
		[factorCovScale{r}, Xr{r}] = Hmsc.scaleMatrix(m.Xr{r}, m.factorCovScaleFlag{r}, m.factorCov(r)-1);
	end
end

if ~isempty(start)
	p = start;
	beta = p.beta;
	gamma = p.gamma;
	iV = inv(p.V);
	if any(m.covScaleFlag~=0)
		ind = m.covScaleFlag==2;
		for i = 1:m.nc
			if m.covScaleFlag(i) == 1
				me = covScale(1,i);
				sd = covScale(2,i);
				if any(ind)
					beta(ind,:) = beta(ind,:) + me*beta(i,:);
					gamma(ind) = gamma(ind) + me*gamma(i);
				end
				beta(i,:) = beta(i,:)*sd;
				gamma(i) = gamma(i)*sd;
				iV(i,:) = iV(i,:)/sd;
				iV(:,i) = iV(:,i)/sd;
			end
		end
	end
	sigma = p.sigma;
	rho = p.rho;
	ph = p.ph;
	nf = p.nf;
	lambda = p.lambda;
	for r=1:m.nr
		if m.factorCov(r)>0 && any(m.factorCovScaleFlag{r}~=0)
			ind = m.factorCovScaleFlag{r}==2;
			for i = 1:m.ncr(r)
				if m.factorCovScaleFlag{r}(i) == 1
					me = factorCovScale{r}(1,i);
					sd = factorCovScale{r}(2,i);
					if any(ind)
						lambda{r}(:,:,ind) = lambda{r}(:,:,ind) + me*lambda{r}(:,:,i);
					end
					lambda{r}(:,:,i) = lambda{r}(:,:,i)*sd;
				end
			end
		end
	end
	eta = p.eta;
	delta = p.delta;
	psijh = p.psijh;
	alpha = p.alpha;
	nfs = p.nfs;
	lambdas = p.lambdas;
	etas = p.etas;
	if m.includeXs > 0 && any(m.sCovScaleFlag~=0)
		ind = m.sCovScaleFlag==2;
		for i = 1:m.ncs
			if m.sCovScaleFlag(i) == 1
				me = sCovScale(1,i);
				sd = sCovScale(2,i);
				if any(ind)
					etas(ind,:) = etas(ind,:) + me*etas(i,:);
				end
				etas(i,:) = etas(i,:)*sd;
			end
		end
	end
	deltas = p.deltas;
	psijhs = p.psijhs;
	%z = Hmsc.update_z(m.Y, X,m.Xr,m.Xs,beta,eta,lambda,etas,lambdas,m.Y,m.pi,m.ncr,m.dist,sigma,m.speciesX,m.includeXs,m.factorCov);
	z = m.Y;
else
	ph = ones([m.ns 1]);
	if m.nr > 0
		% 		nf = ones(m.nr,1)*max(2,floor(log(m.ns)*3));        % initial number of factors. k(1) for residual, k(i1+1) for random effects
		nf=ones(m.nr,1);
		if m.ns>1
			nf = nf.*2;
		end
	else
		nf = [];
	end
	if m.ncs > 0
		nfs = 2;
	else
		nfs = 0;
	end
	sigma = eye(m.ns);
	beta = zeros(m.nc,m.ns);
	warning('off','all');
	with_data = 1-isnan(m.Y);
	for j = 1:m.ns
		if m.speciesX
			X1 = X{j};
		else
			X1 = X;
		end
		if sum(with_data(:,j))>m.nc+1
			if or(m.dist(j,1)==1, m.dist(j,1)==4)
				beta(:,j) = glmfit(X1,m.Y(:,j),'normal','constant','off');
				if m.dist(j,2)
					res=X1*beta(:,j)-m.Y(:,j);
					sigma(j,j)=max(0.0001,(res'*res)/m.ny);
				elseif m.dist(j,1)==4
					sigma(j,j)=0.01;
				end
			end
			if m.dist(j,1)==2
				beta(:,j) = glmfit(X1,m.Y(:,j),'binomial','link','probit','constant','off');			
			end
			if m.dist(j,1)==3
				beta(:,j) = glmfit(X1,log(m.Y(:,j)+0.5),'normal','constant','off');
				if m.dist(j,2)
					res=X1*beta(:,j)-log(m.Y(:,j)+0.5);
					sigma(j,j)=max(0.0001,(res'*res)/m.ny);
				else
					sigma(j,j)=0.01;
				end
			end
		end
	end
	beta=max(min(beta,3),-3);
	warning('on','all');
	Vn = inv(m.T'*m.T +eye(m.nt));
	gamma = Vn*(m.T'*beta');
	be0 = gamma'*m.T';
	V = cov(beta')+0.1*eye(m.nc);
	iV = inv(V);
	
	rho=1;
	
	if m.nr==0
		alpha=[];
		lambda=[];
		eta=[];
		delta=[];
		psijh=[];
	end
	eta = cell(1, m.nr);
	alpha = cell(1, m.nr);
	for i=1:m.nr
		eta1 = normrnd(0,1,[m.np(i),nf(i)]);
		eta{i}=eta1;
		if m.spatial(i)
			alpha1 = round(size(m.alphapw{i},1)/2).*ones(nf(i),1);
			alpha{i}=alpha1;
		else
			alpha{i} = [];
		end
	end
	
	etas=normrnd(0,1,[m.ncs,nfs]);
	lambdas=normrnd(0,1,[nfs,m.ns]);
	
	for i1=1:m.nr
		if m.factorCov(i1)
			psijh1 = gamrnd(m.nur/2,2/m.nur,[m.ns,nf(i1),m.ncr(i1)]);
			delta1 = [gamrnd(m.a1r,1/m.b1r,[1,m.ncr(i1)]);gamrnd(m.a2r,1/m.b2r,[nf(i1)-1,m.ncr(i1)])];
			lambda1 = normrnd(0,1,[nf(i1),m.ns,m.ncr(i1)]);
		else
			psijh1 = gamrnd(m.nur/2,2/m.nur,[m.ns,nf(i1)]);
			delta1 = [gamrnd(m.a1r,1/m.b1r);gamrnd(m.a2r,1/m.b2r,[nf(i1)-1,1])];
			lambda1 = normrnd(0,1,[nf(i1),m.ns]);
		end
		psijh{i1} = psijh1;
		delta{i1} = delta1;
		lambda{i1}=lambda1;
	end
	
	z =  Hmsc.update_z(m.Y,X,m.Xr,m.Xs,beta,eta,lambda,etas,lambdas,m.Y,m.pi,m.ncr,m.dist,sigma,m.speciesX,m.includeXs,m.factorCov);
	z(isnan(z))=0;
end

psijhs = gamrnd(m.nus/2,2/m.nus,[m.ns,nfs]);
deltas = [gamrnd(m.a1s,1/m.b1s);gamrnd(m.a2s,1/m.b2s,[nfs-1,1])];
end


