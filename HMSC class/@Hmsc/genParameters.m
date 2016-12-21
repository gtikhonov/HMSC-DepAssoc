function genParameters(m)

% add console output Generating XXXXX

if isempty(m.truePar)
	p = HmscPar();
else
	p = m.truePar;
end;

if isempty(p.nf)
	p.nf = 2*ones(1, m.nr);
end

if isempty(p.sigma)
	p.sigma = eye(m.ns);
	for j = 1:m.ns
		if (m.dist(j,1)==3 || m.dist(j,1)==4) && ~m.dist(j,2)
			p.sigma(j,j)=0.01;
		elseif m.dist(j,1)~= 2
			p.sigma(j,j) = 1/gamrnd(m.asigma(j),1/m.bsigma(j) );
		end
	end
end

if isempty(p.gamma)
	p.gamma = normrnd(0, 1, [m.nt, m.nc]);
	p.gamma(1) = 0;
end

if isempty(p.V)
	df = m.nc+2;
	p.V = iwishrnd(m.V0,df)*(df-m.nc-1);
end;

if isempty(p.ph)
	p.ph = ones(m.ns, 1);
	if m.outlierSpecies
		for j = 1:min(5,round(m.ns/10))
			p.ph(j) = 0.0001;
		end
	end
end

if m.phylogeny
	if isempty(p.rho)
		rhop = m.rhopw(:,1);
		rhow = m.rhopw(:,2);
		p.rho = randsample(1:length(rhop), 1, true, rhow);
	end
end

if isempty(p.beta)
	mbeta = m.T*p.gamma;
	if m.phylogeny
		rhop = m.rhopw(:,1);
		if rhop(p.rho) >= 0
			rhoC = rhop(p.rho)*m.C;
		else
			rhoC = (-rhop(p.rho))*corrcov(inv(m.C));
		end
		Q = rhoC+(1-abs(rhop(p.rho)))*eye(m.ns);
		M = kron(p.V,Q);
		if m.outlierSpecies
			PHI = (1./sqrt(p.ph))*(1./sqrt(p.ph))';
			MM = kron(ones(m.nc),PHI).*M;
		else
			MM = M;
		end
		p.beta = mvnrnd(mbeta(:), MM)';
		p.beta = reshape(p.beta, m.ns, m.nc)';
	else
		p.beta = nan(m.ns, m.nc);
		for j = 1:m.ns
			p.beta(j,:) = mvnrnd(mbeta(j,:), (1/p.ph(j))*p.V);
		end
		p.beta = p.beta';
	end
end

if isempty(p.alpha) % all( cellfun( @isempty, p.alpha) )
	p.alpha = cell(1, m.nr);
	for i = 1:m.nr
		alpha1 = NaN;
		if m.spatial(i)
			xy1 = m.xy{i};
			alphap = m.alphapw{i}(:,1);
			alphaw = m.alphapw{i}(:,2);
			alpha1 = randsample(1:length(alphap), p.nf(i), true, alphaw);
		end
		p.alpha{i} = alpha1;
	end
end

if isempty(p.eta) %all( cellfun( @isempty, p.eta) )
	p.eta = cell(1, m.nr);
	for i = 1:m.nr
		if m.spatial(i)
			alpha1 = p.alpha{i};
			alphap = m.alphapw{i}(:,1);
			xy1 = m.xy{i};
			di = zeros(m.np(i));
			for j = 1:m.spatDim(i)
				xx = repmat(xy1(:,j), 1, m.np(i));
				dx = xx-xx';
				di = di+dx.^2;
			end
			di = sqrt(di);
			eta1=[];
			for h=1:p.nf(i)
				if(alphap(alpha1(h))>10^(-7))
					W = exp(-di/alphap(alpha1(h)));
				else
					W = eye(length(di));
				end
				eta1 = [eta1; mvnrnd(zeros(m.np(i),1),W)];
			end
			eta1 = eta1';
		else
			eta1 = normrnd(0, 1, [m.np(i), p.nf(i)]);
		end
		p.eta{i} = eta1;
	end
end

if isempty(p.delta) %isempty( p.delta )
	p.delta = cell(1, m.nr);
	for i = 1:m.nr
		if m.factorCov(i)
			delta1 = zeros(p.nf(i), m.ncr(i));
			delta1(1,:) = gamrnd(m.a1r, 1/m.b1r, [1,m.ncr(i)]);
			for h = 2:p.nf(i)
				delta1(h,:) = gamrnd(m.a2r, 1/m.b2r, [1,m.ncr(i)]);
			end
		else
			delta1=zeros(p.nf(i),1);
			delta1(1) = gamrnd(m.a1r,1/m.b1r);
			for h = 2:p.nf(i)
				delta1(h) = gamrnd(m.a2r,1/m.b2r);
			end
		end
		p.delta{i} = delta1;
	end
end

if isempty(p.psijh) %all( cellfun( @isempty, p.psijh ) )
	p.psijh = cell(1, m.nr);
	for i = 1:m.nr
		if m.factorCov(i)
			psijh1 = gamrnd(m.nur/2, 1./(m.nur/2), [p.nf(i),m.ns,m.ncr(i)]);
		else
			psijh1 = gamrnd(m.nur/2, 1./(m.nur/2), p.nf(i), m.ns);
		end
		p.psijh{i} = permute(psijh1,[2,1,3]);
	end
end

if isempty(p.lambda) %all( cellfun( @isempty, p.lambda ) )
	p.lambda = cell(1, m.nr);
	for i = 1:m.nr
		delta1 = p.delta{i};
		psijh1 = permute(p.psijh{i},[2,1,3]);
		tauh1 = cumprod(delta1);
		if m.factorCov(i)
			tauh2 = permute(tauh1, [1,3,2]);
			lambda1 = normrnd(0, sqrt(1./(repmat(tauh2,1,m.ns,1).*psijh1)));
		else
			lambda1=normrnd(0, sqrt(1./(repmat(tauh1,1,m.ns).*psijh1)));
		end
		p.lambda{i} = lambda1;
	end
end


% for i = 1:m.nr
% 	if m.factorCov(i)
% 		delta1 = zeros(p.nf(i), m.ncr(i));
% 		delta1(1,:) = gamrnd(m.a1r, 1/m.b1r, [1,m.ncr(i)]);
% 		for h = 2:p.nf(i)
% 			delta1(h,:) = gamrnd(m.a2r, 1/m.b2r, [1,m.ncr(i)]);
% 		end
% 		tauh1 = cumprod(delta1);
% 		tauh2 = permute(tauh1, [1,3,2]);
% 		psijh1 = gamrnd(m.nur/2, 1./(m.nur/2), [p.nf(i),m.ns,m.ncr(i)]);
% 		lambda1 = normrnd(0, sqrt(1./(repmat(tauh2,1,m.ns,1).*psijh1)));
% 	else
% 		delta1=zeros(p.nf(i),1);
% 		delta1(1) = gamrnd(m.a1r,1/m.b1r);
% 		for h = 2:p.nf(i)
% 			delta1(h) = gamrnd(m.a2r,1/m.b2r);
% 		end
% 		tauh1 = cumprod(delta1);
% 		psijh1 = gamrnd(m.nur/2, 1./(m.nur/2), p.nf(i), m.ns);
% 		lambda1=normrnd(0, sqrt(1./(repmat(tauh1,1,m.ns).*psijh1)));
% 	end
% 	p.lambda{i} = lambda1;
%
% 	p.delta{i} = delta1;
% 	p.psijh{i} = permute(psijh1,[2,1,3]);
% end

if m.includeXs
	if isempty(p.nfs)
		p.nfs = 2;
	end
	if isempty(p.etas)
		p.etas = normrnd(0,1,[m.ncs,p.nfs]);
	end
	if isempty(p.lambdas)
		p.lambdas = normrnd(0,1,[p.nfs,m.ns]);
	end
end

m.truePar = p;

end