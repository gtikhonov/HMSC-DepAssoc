function p = addToOutput(m,beta,gamma,sigma,rho,iV,lambda,eta,alpha,ph,nf,lambdas,etas,nfs,psijh,delta,psijhs,deltas,covScale,sCovScale,factorCovScale)
p = HmscPar();
if any(m.covScaleFlag~=0)
	ind = m.covScaleFlag==2;
	for i = 1:m.nc
		if m.covScaleFlag(i) == 1
			me = covScale(1,i);
			sd = covScale(2,i);
			beta(i,:) = beta(i,:)/sd;
			if any(ind)
				beta(ind,:) = beta(ind,:) - me*beta(i,:);
			end
			gamma(i) = gamma(i)/sd;
			iV(i,:) = iV(i,:)*sd;
			iV(:,i) = iV(:,i)*sd;
		end
	end
end
if m.includeXs > 0 && any(m.sCovScaleFlag~=0)
	ind = m.sCovScaleFlag==2;
	for i = 1:m.ncs
		if m.sCovScaleFlag(i) == 1
			me = sCovScale(1,i);
			sd = sCovScale(2,i);
			etas(i,:) = etas(i,:)/sd;
			if any(ind)
				etas(ind,:) = etas(ind,:) - me*etas(i,:);
			end
		end
	end
end
for r=1:m.nr
	if m.factorCov(r)>0 && any(m.factorCovScaleFlag{r}~=0)
		ind = m.factorCovScaleFlag{r}==2;
		for i = 1:m.ncr(r)
			if m.factorCovScaleFlag{r}(i) == 1
				me = factorCovScale{r}(1,i);
				sd = factorCovScale{r}(2,i);
				lambda{r}(:,:,i) = lambda{r}(:,:,i)/sd;
				if any(ind)
					lambda{r}(:,:,ind) = lambda{r}(:,:,ind) - me*lambda{r}(:,:,i);
				end
			end
		end
	end
end

p.beta = beta;
p.gamma = gamma;
p.sigma = sigma;
p.V = inv(iV);
p.lambda = lambda;
p.eta = eta;
p.ph = ph;
p.nf = nf;
p.lambdas = lambdas;
p.nfs = nfs;
p.etas = etas;

p.psijh = psijh;
p.delta = delta;
p.psijhs = psijhs;
p.deltas = deltas;


p.rho = rho;
p.alpha = alpha;

end