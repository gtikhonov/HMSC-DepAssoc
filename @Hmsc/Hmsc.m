classdef Hmsc < handle
	properties
		folder
	end
	properties (SetAccess = private)
		% all properties
		% model definitions - setted with constructor
		traits,speciesX,phylogeny,includeXs,outlierSpecies, speciesXs
		spatial,factorCov
		% obligatory dimensions
		ny, ns, nc
		% optional dimensions
		nt, ncr, ncs, spatDim
		% obligatory data
		X, dist, Y, pi, piCell
		% optional data
		T, Xr, XrCell, Xs, xy, xyCell, C
		% defined from other data data
		nr, np
		iWg, detWg
		spNames, covNames, traitNames, sCovNames, levelNames
		covScale, covScaleFlag, sCovScale, sCovScaleFlag, factorCovScale, factorCovScaleFlag;
		
		% priors
		nu
		nur, a1r, a2r, b1r, b2r
		nus, a1s, a2s, b1s, b2s
		asigma, bsigma
		f0, V0, Ugamma, mgamma
		alphapw
		rhopw
		
		% true parameters - paremeters are incorporated to a separate simple class
		truePar
		% rep parameters - cell array of parameters replicates
		repPar, repN
		
		% sampler options
		samples, thinning,
		adapt, fixNf, adaptXs, fixNfs
		samVec
		start
		ramPost, stfPost
		
		% post posterior sampling options
		postSamVec, postSamInd, postSamN
		
		piKey, piMap;
	end
	properties (Access = private)
	end
	methods
		function m = Hmsc(folder, traits, speciesX, phylogeny, includeXs, outlierSpecies, spatial, factorCov)
			if length(spatial) ~= length(factorCov)
				error('HMSC: vectors spatial and factorCov must be of equal length');
			end
			m.folder = folder;
			m.traits = traits;
			m.speciesX = speciesX;
			m.phylogeny = phylogeny;
			m.includeXs = includeXs;
			m.outlierSpecies = outlierSpecies;
			m.spatial = spatial;
			m.factorCov = factorCov;
			m.nr = length(spatial);
			m.spatDim = NaN(1, m.nr);
			if any(m.factorCov) == false
				m.ncr = ones(1, m.nr);
			end
			if m.includeXs == 0
				m.ncs = 0;
			end
		end
		
		%Set model data
		setData(m, Y, dist, X, piCell, xyCell, T, C, XrCell, Xs)
		setDim(m, ny, ns, nc)
		setDist(m, dist)
		setX(m, X)
		setPi(m, piCell)
		setY(m, Y)
		setTraits(m, T)
		setSpatialLocations(m, xyCell)
		setRanFactorCov(m, XrCell)
		setPhylogeny(m, C)
		setXs(m, Xs)
		setFolder(m, folder)
		
		setPriorsDefault(m)
		setPriorsRandomFactors(m, nur, alr, b1r, a2r, b2r)
		setPriorsXs(m, nus, als, b1s, a2s, b2s)
		setPriorsSigma(m, asigma, bsigma)
		setPriorsGamma(m, f0, V0, Ugamma, mgamma)
		setPriorsOutlierSpecies(m, nu)
		setPriorsRho(m, rhop, rhow)
		setPriorsAlpha(m, alphap, alphaw)
		setParameters(m, par)
		
		%set names
		setSpeciesNames(m, spNames)
		setCovNames(m, covNames)
		setSCovNames(m, sCovNames)
		setTraitNames(m, traitNames)
		setLevelNames(m, levelNames)
		
		setCovScaling(m, scaleFlagVec)
		
		% data radomizers
		genX(m) 
		genTraits(m, nt) 
		genPhylogeny(m) 
		genPi(m, np) 
		genRanFactorCov(m, ncr)
		genSpatialLocations(m, spatDim) 
		genXs(m, ncs)
		genY(m, misFraction) 
		genParameters(m) 
		
		% sampler
		setMCMCOptions(m, samples, thinning)
		setMCMCAdapt(m, adapt, fixNf, adaptXs, fixNfs)
		setMCMCSaveOptions(m, ramPost, stfPost)
		sampleMCMC(m, nRuns, append, startPar, verb)
		
		postFileToRam(m, loadVec)
		postRamToFile(m, loadVec)
		postRamClear(m)
		postFileClear(m)
		postRemove(m)
		
		createPostSamVec(m, runs, samVec);
		setPostThinning(m, runs, postThin)
		
		summary(m, par, quantiles, dispTrue, showToScreen, saveToFile);
		
		plotGamma(m, dispTrue, showToScreen, saveToFile, type)
		plotBeta(m, dispTrue, showToScreen, saveToFile, type)
		plotV(m, dispTrue, showToScreen, saveToFile, type)
		plotSigma(m, dispTrue, showToScreen, saveToFile, type)
		plotRho(m, dispTrue, showToScreen, saveToFile, type)
		plotPh(m, dispTrue, showToScreen, saveToFile, type)
		plotLambda(m, dispTrue, showToScreen, saveToFile, type)
		plotEta(m, dispTrue, showToScreen, saveToFile, type)
		plotAlpha(m, dispTrue, showToScreen, saveToFile, type)
		plotOmega(m, Xr, maxPairs, dispTrue, showToScreen, saveToFile, type)
		plotEtas(m, dispTrue, showToScreen, saveToFile, type)
		plotLambdas(m, dispTrue, showToScreen, saveToFile, type)
		plotAm(m, dispTrue, showToScreen, saveToFile, type)
		plotAll(m, dispTrue, type)
		
		res = getPostGamma(m)
		res = getPostBeta(m)
		res = getPostV(m)
		res = getPostSigma(m)
		res = getPostRho(m)
		res = getPostLambda(m, i, k)
		res = getPostOmega(m, i, X)
		res = getPostEta(m, i)
		res = getPostAlpha(m, i)
		res = getPostAm(m)
		res = getPostEtas(m)
		res = getPostLambdas(m)
		res = getPostPh(m)
		
		Y = predict(m, n, X, piCell, xyCell, XrCell, expected) % X and pi must be specified, everything else - only if applicable.
		Y = predictConditional(m, n, Yc, nmcmc, X, piCell, xyCell, XrCell, expected) % X and pi must be specified, everything else - only if applicable.
		
		[correlations,support,index] = computeCorrelations(m,level,Xr,threshold)
		R2 = computeR2(m,predN)
		[fixed, fixedsplit, random, traitR2] = computeVariances(m,group);
		plotR2(m,R2,showToScreen, saveToFile)
		plotVariancePartitioning(m,fixed,fixedsplit,random,traitR2,groupnames,grouprandomnames,showToScreen, saveToFile)
		plotCorrelations(m,correlations, index, plottitle, type, showToScreen, saveToFile)
	end
	
	methods (Access = private)
		[iQg,detQg,rhopw,iWgA,detWgA,alphapwA] = calculateDataParameters(m);
		[X,covScale,Xs,sCovScale,Xr,factorCovScale,beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs] = computeInitialValues(m, start);
		stfMCMCRun(m, parVec, run)
		summarize(m, name, values, valuesT, q, d1, d2, tran, rNames, cNames, dispTrue, showToScreen, saveToFile);
		p = addToOutput(m,beta,gamma,sigma,rho,iV,lambda,eta,alpha,ph,nf,lambdas,etas,nfs,psijh,delta,psijhs,deltas,covScale,sCovScale,factorCovScale)
		
		setMCMCVector(m, samVec);
		
		plotMixing(m, values,valuesT,dispTrue,label, showToScreen, saveToFile)
		plotBox(m, values, valuesT, dispTrue, label, showToScreen, saveToFile)
	end
	
	methods (Static = true)
		lambda = update_lambda(X,Xr,Xs,z,beta,sigma,eta,lambda,etas,lambdas,delta,psijh,pi,nf,ncr,spatial,factorCov,speciesX,includeXs)
		eta = update_eta(z,X,Xr,Xs,beta,sigma,eta,alpha,lambda,etas,lambdas,nf,pi,ncr,spatial,factorCov,iWg,speciesX,includeXs)
		alpha = update_alpha(eta,alpha,nf,pi,spatial,iWg,detWg,alphapw)
		z = update_z(z, X,Xr,Xs,beta,eta,lambda,etas,lambdas,Y,pi,ncr,dist,sigma,speciesX,includeXs,factorCov)
		sigma = update_sigma(X,Xr,Xs,sigma,beta,eta,lambda,etas,lambdas,pi,ncr,dist,asigma,bsigma,z,spatial,factorCov,speciesX,includeXs)
		beta = update_beta(X,Xr,Xs,z,pi,ncr,eta,lambda,etas,lambdas,gamma,sigma,T,ph,iV,rho,phylogeny,iQg,detQg,outlierspecies,spatial,factorCov,speciesX,includeXs)
		[gamma,iV] = update_gamma_V(T,beta,gamma,ph,rho,V0,f0,Ugamma,mgamma,phylogeny,iQg,outlierspecies)
		rho = update_rho(rho,beta,T,gamma,detQg,iQg,iV,rhopw,ph,outlierspecies)
		[psijh,delta] = update_lambda_priors(nf,nur,a1r,a2r,b1r,b2r,psijh,delta,lambda,factorCov)
		etas=update_etas(X,Xs,Xr,pi,ncr,z,lambdas,beta,sigma,eta,lambda,spatial,speciesX,includeXs,factorCov)
		lambdas=update_lambdas(X,Xs,Xr,pi,ncr,z,lambdas,etas,beta,sigma,eta,lambda,psijhs,deltas,spatial,speciesX,includeXs,factorCov)
		[psijhs,deltas] = update_lambdas_priors(nfs,nus,a1s,a2s,b1s,b2s,psijhs,deltas,lambdas)
		[nf,lambda,eta,alpha,psijh,delta] = update_nf(nf,ncr,ns,mcmc,np,nur,a2r,b2r,lambda,eta,alpha,psijh,delta,spatial)
		[nfs,lambdas,etas,psijhs,deltas] = update_nfs(mcmc,nus,a2s,b2s,lambdas,etas,psijhs,deltas)
		ph = update_ph(T,beta,iV,ph,gamma,nu)

		eps = tnormrnd(mu,si,low,high)		
		[scale, XSc] = scaleMatrix(X, scaleFlag, speciesX)
	end
end
