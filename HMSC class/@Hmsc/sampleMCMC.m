function sampleMCMC(m, nRuns, append, startPar, verb)

speciesX = m.speciesX;
phylogeny = m.phylogeny;
includeXs = m.includeXs;
outlierSpecies = m.outlierSpecies;
spatial = m.spatial;
factorCov = m.factorCov;

f0 = m.f0;
V0 = m.V0;
Ugamma = m.Ugamma;
mgamma = m.mgamma;
asigma = m.asigma;
bsigma = m.bsigma;
nur = m.nur;
a1r = m.a1r;
b1r = m.b1r;
a2r = m.a2r;
b2r = m.b2r;
nus = m.nus;
a1s = m.a1s;
b1s = m.b1s;
a2s = m.a2s;
b2s = m.b2s;
nu = m.nu;

% nc = m.nc;
ns = m.ns;
% nt = m.nt;
ncr = m.ncr;
np = m.np;
dist = m.dist;

Xr = m.Xr;
Xs = m.Xs;
T = m.T;
Y = m.Y;
pi = m.pi;
rhopw = m.rhopw;
alphapw = m.alphapw;

[iQg,detQg,iWg,iWgR,detWg] = m.calculateDataParameters();
m.iWg = iWg;
m.detWg = detWg;

if append == false
	m.repPar = [];
	m.repN = 0;
	startRunN = 0;
	[X,covScale,Xs,sCovScale,Xr,factorCovScale, beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs] = m.computeInitialValues(startPar);
	% clear the directory
else
	startRunN = length(m.repPar);
	flag = 0;
	if isempty( m.repPar{m.repN} )
		flag = 1;
		m.postFileToRam(m.repN);
	end
	lastPar = m.repPar{m.repN}{end};
	if flag
		m.postRamClear();
	end
	[X,covScale,Xs,sCovScale,Xr,factorCovScale,beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs] = m.computeInitialValues(lastPar);
end
tic;
if verb > 0
	fprintf('Sampling starts.\n');
end
for run = startRunN+(1:nRuns)
	if verb > 1
		fprintf(strcat('SAVE REPLICATE:',int2str(run),'\n'));
	end
	parVec = cell(1, length(m.samVec) );
	parK = 1;
	for mcmc = 1:m.samVec(end)
		lambda = Hmsc.update_lambda(X,Xr,Xs,z,beta,sigma,eta,lambda,etas,lambdas,delta,psijh,pi,nf,ncr,spatial,factorCov,speciesX,includeXs);
		eta = Hmsc.update_eta(z,X,Xr,Xs,beta,sigma,eta,alpha,lambda,etas,lambdas,nf,pi,ncr,spatial,factorCov,iWg,speciesX,includeXs);
		alpha = Hmsc.update_alpha(eta,alpha,nf,pi,spatial,iWg,iWgR,detWg,alphapw);
		z = Hmsc.update_z(z, X,Xr,Xs,beta,eta,lambda,etas,lambdas,Y,nf,pi,ncr,dist,sigma,speciesX,includeXs,factorCov);
		sigma = Hmsc.update_sigma(X,Xr,Xs,sigma,beta,eta,lambda,etas,lambdas,nf,pi,ncr,dist,asigma,bsigma,z,spatial,factorCov,speciesX,includeXs);
		beta = Hmsc.update_beta(X,Xr,Xs,z,nf,pi,ncr,eta,lambda,etas,lambdas,gamma,sigma,T,ph,iV,rho,phylogeny,iQg,detQg,outlierSpecies,spatial,factorCov,speciesX,includeXs);
		[gamma,iV] = Hmsc.update_gamma_V(T,beta,gamma,ph,rho,V0,f0,Ugamma,mgamma,phylogeny,iQg,outlierSpecies);
		if phylogeny
			rho = Hmsc.update_rho(rho,beta,T,gamma,detQg,iQg,iV,rhopw,ph,outlierSpecies);
		end
		if outlierSpecies
			ph = Hmsc.update_ph(T,beta,iV,ph,gamma,nu);
		end
		[psijh,delta] = Hmsc.update_lambda_priors(nf,nur,a1r,a2r,b1r,b2r,psijh,delta,lambda,factorCov);
		
		if includeXs > 0
			etas = Hmsc.update_etas(X,Xs,Xr,pi,ncr,z,lambdas,beta,sigma,eta,lambda,spatial,speciesX,includeXs,factorCov);
			lambdas = Hmsc.update_lambdas(X,Xs,Xr,pi,ncr,z,lambdas,etas,beta,sigma,eta,lambda,psijhs,deltas,spatial,speciesX,includeXs,factorCov);
			[psijhs,deltas] = Hmsc.update_lambdas_priors(nfs,nus,a1s,a2s,b1s,b2s,psijhs,deltas,lambdas);
			if (run<m.adaptXs(1) || (run==m.adaptXs(1) && mcmc<=m.adaptXs(2)) ) && m.fixNfs==false
				[nfs,lambdas,etas,psijhs,deltas] = Hmsc.update_nfs(mcmc,nus,a2s,b2s,lambdas,etas,psijhs,deltas);
			end
		end
		
		if m.nr > 0
			if (run<m.adapt(1) || (run==m.adapt(1) && mcmc<=m.adapt(2)) ) && isequal(m.fixNf, false)
				[nf,lambda,eta,alpha,psijh,delta] = Hmsc.update_nf(nf,ncr,ns,mcmc,np,nur,a2r,b2r,lambda,eta,alpha,psijh,delta,spatial);
			end
		end
		
		if any(mcmc==m.samVec)
			if verb > 2
				fprintf(strcat('save replicate:', int2str(run), ',run:', int2str(run-startRunN), ',iteration round:',int2str(mcmc),'---', num2str(round(toc)),'\n'));
			end
			p = m.addToOutput(beta,gamma,sigma,rho,iV,lambda,eta,alpha,ph,nf,lambdas,etas,nfs,psijh,delta,psijhs,deltas,covScale,sCovScale,factorCovScale);
			parVec{parK} = p;
			parK = parK + 1;
		end
	end
	if m.ramPost == true
		m.repPar{run} = parVec;
	else
		m.repPar{run} = [];
	end
	if m.stfPost == true
		m.stfMCMCRun(parVec, run);
	end
	m.repN = m.repN + 1;
end
if verb>0
	toc;
end
end