function Y = predictConditional(m, n, Yc, nmcmc, X, piCell, xyCell, XrCell, expected)

Xs = m.Xs;

piN = nan(size(piCell));
piKeyN = cell(1, m.nr);
piMapN = cell(1, m.nr);
xy = cell(1, m.nr);
Xr = cell(1, m.nr);
for r=1:m.nr
	piKey = m.piKey{r};
	keySet = unique(piCell(:,r));
	exInd = m.piMap{r}.isKey(keySet);
	if any(~exInd)
		keySetN = keySet(~exInd);
		piKeyN{r} = [piKey; keySetN];
		newMap = containers.Map(keySetN, (length(piKey)+1):length(piKeyN{r}) );
		piMapN{r} = [m.piMap{r}; newMap];
	else
		piKeyN{r} = piKey;
		piMapN{r} =  m.piMap{r};
	end
	piN(:,r) = cell2mat(piMapN{r}.values( piCell(:,r)) );
	if m.spatial(r)
		piKey = piKeyN{r};
		xyKey = xyCell{r}(:,1);
		xyMap = containers.Map(xyKey,1:length(xyKey));
		if ~all(xyMap.isKey(piKey))
			error('HMSC: some units, defined at level %d were not given spatial coordinates\n', r);
		end
		ind = cell2mat( xyMap.values(piKey) );
		xy1 = xyCell{r}(ind, 2:size(xyCell{r}, 2));
		xy1 = cell2mat(xy1);
		xy{r} = xy1;
	end
	if m.factorCov(r)
		piKey = piKeyN{r};
		XrKey = XrCell{r}(:,1);
		XrMap = containers.Map(XrKey,1:length(XrKey));
		if ~all(XrMap.isKey(piKey))
			error('HMSC: some units, defined at level %d were not given spatial coordinates\n', r);
		end
		ind = cell2mat( XrMap.values(piKey) );
		Xr1 = XrCell{r}(ind, 2:size(XrCell{r}, 2));
		Xr1 = cell2mat(Xr1);
		Xr{r} = Xr1;
	end
end
pi = piN;

res = cell(1, n);
r = randi([1, m.postSamN],1,n);
dis = cell(1, m.nr);
for i = 1:m.nr
	newPiInd = ~ismember(pi(:,i), m.pi(:,i));
	pi1 = pi(newPiInd, i);
	newPiN = length( unique(pi1) );
	di = [];
	if m.spatial(i) && newPiN > 0
		di = zeros(m.np(i)+newPiN);
		xy1 = xy{i};
		for j = 1:m.spatDim(i)
			xx = repmat(xy1(:,j), 1, m.np(i)+newPiN);
			dx = xx-xx';
			di = di+dx.^2;
		end
	end
	di = sqrt(di);
	dis{i} = di;
end

for rN = 1:n
	if ~mod(rN, 10)
		fprintf('Calculating conditional prediction %d\n', rN);
	end
	p = m.postSamVec(r(rN));
	if m.speciesX
		ny = size(X{1}, 1);
		Ez = zeros(ny, m.ns);
		for i = 1:m.ns
			Ez(:,i) = X{i}*p.beta(:,i);
		end
	else
		ny = size(X, 1);
		Ez = X*p.beta;
	end
	
	eta = cell(1, m.nr);
	for i = 1:m.nr
		etaM = p.eta{i};
		newPiInd = ~ismember(pi(:,i), m.pi(:,i));
		pi1 = pi(newPiInd, i);
		newPiN = length( unique(pi1) );
		if m.spatial(i) && newPiN > 0
			di = dis{i};
			alphaInd = p.alpha{i};
			alphapw = m.alphapw{i};
			etaN = zeros(newPiN, p.nf(i));
			for j = 1:p.nf(i)
				alpha = alphapw(alphaInd(j), 1);
				if alpha > 0
					% 					di11 = di(1:m.np(i), 1:m.np(i));
					% 					W11 = exp(-di11/alpha);
					di12 = di(1:m.np(i), m.np(i)+(1:newPiN));
					di22 = di(m.np(i)+(1:newPiN), m.np(i)+(1:newPiN));
					W12 = exp(-di12/alpha);
					W22 = exp(-di22/alpha);
					iW11 = m.iWg{i}(:,:,alphaInd(j));
					muN = W12' * iW11 * etaM(:,j);
					WN = W22 - W12' * iW11 * W12;
					WN = (WN+WN') / 2;
					etaN(:, j) = mvnrnd(muN, WN);
				else
					etaN(:, j) = normrnd(0,1,[newPiN,1]);
				end
			end
		else
			etaN = normrnd(0,1,[newPiN,p.nf(i)]);
		end
		eta{i} = [etaM; etaN];
	end
	
	z=zeros(ny,m.ns);
	if nmcmc>1
		for level=1:m.nr
			iWg=[];
			nalphai=[];
			if m.spatial(level)
				di=dis{level};
				nalphai=1:p.nf(level);
				for ag=1:p.nf(level)
					alp = m.alphapw{level}(p.alpha{level}(ag),1);
					if alp < 1e-5
						W = eye(size(di,1));
					else
						W = exp(-di/alp);
					end
					iWg = cat( 3, iWg, inv(W) );
				end
			end
			nalpha{level}=nalphai;
			iWgA{level}=iWg;
		end
	end
	for i=1:nmcmc
		z=m.update_z(z,X,Xr,Xs,p.beta,eta,p.lambda,p.etas,p.lambdas,Yc,pi,m.ncr,m.dist,p.sigma,m.speciesX,m.includeXs,m.factorCov);
		if (i<nmcmc)
			eta=m.update_eta(z,X,Xr,Xs,p.beta,p.sigma,eta,nalpha,p.lambda,p.etas,p.lambdas,p.nf,pi,m.ncr,m.spatial,m.factorCov,iWgA,m.speciesX,m.includeXs);
		end
	end
	
	Y = z;
	for j = 1:m.ns
		if(m.dist(j,1) == 2)
			if expected
				Y(:,j) = normcdf(z(:,j));
			else
				Y(:,j) = z(:,j)>0;
			end
		end
		if(m.dist(j,1) == 3)
			if expected
				Y(:,j) = exp(z(:,j));
			else
				Y(:,j) = floor(exp(z(:,j)));
			end
		end
		if(m.dist(j,1) == 4)
			if expected
				Y(:,j) = max(0,z(:,j));
			else
				Y(:,j) = max(0,floor(z(:,j)));
			end
			Y(:,j) = max(0,floor(z(:,j)));
		end
	end
	res{rN} = Y;
end
Y = res;
end
