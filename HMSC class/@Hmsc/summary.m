function [resMean, varargout] = summary(m, par, quantiles, dispTrue, showToScreen, saveToFile)

if dispTrue && isempty(m.truePar)
	error('HMSC: true parameter set is empty, but was ordered to be shown');
end

if nargout > 0
	if length(quantiles)+dispTrue < nargout-1
		error('HMSC: too many output arguments');
	elseif length(quantiles)+dispTrue > nargout-1
		error('HMSC: too few output arguments');
	end
end



if strcmp(par, 'beta')
	res = m.getPostBeta();
	resCell = m.summarize('beta', res{1}, res{2}, quantiles, m.nc, m.ns, true, m.spNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'gamma')
	res = m.getPostGamma();
	resCell = m.summarize('gamma', res{1}, res{2}, quantiles, m.nt, m.nc, false, m.traitNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'V')
	res = m.getPostV();
	resCell = m.summarize('V', res{1}, res{2}, quantiles, m.nc, m.nc, false, m.covNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'sigma')
	res = m.getPostSigma();
	resCell = m.summarize('sigma', res{1}, res{2}, quantiles, m.ns, 1, false, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'rho')
	res = m.getPostRho();
	resCell = m.summarize('rho', res{1}, res{2}, quantiles, 1, 1, false, {'rho'}, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'ph')
	res = m.getPostPh();
	resCell = m.summarize('ph', res{1}, res{2}, quantiles, m.ns, 1, false, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'etas')
	res = m.getPostEtas();
	resCell = m.summarize('etas', res{1}, res{2}, quantiles, m.ncs, m.postSamVec(1).nfs, false, m.sCovNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'lambdas')
	res = m.getPostLambdas();
	resCell = m.summarize('lambdas', res{1}, res{2}, quantiles, m.postSamVec(1).nfs, m.ns, true, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'Am')
	res = m.getPostAm();
	resCell = m.summarize('Am', res{1}, res{2}, quantiles, m.ncs, m.ns, true, m.spNames, m.sCovNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'alpha')
	resCell = cell(m.nr, 1+length(quantiles)+dispTrue);
	for r=1:m.nr
		if m.spatial(r)
			res = m.getPostAlpha(r);
			label = strcat( 'alpha_', int2str(r) );
			resCell(r,:) = m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(r), 1, false, [], [], dispTrue, showToScreen, saveToFile);
		end
	end
elseif strcmp(par, 'eta')
	resCell = cell(m.nr, 1+length(quantiles)+dispTrue);
	for k = 1:m.nr
		res = m.getPostEta(k);
		label = strcat('eta_', int2str(k));
		resCell(r,:) = m.summarize(label, res{1}, res{2}, quantiles, m.np(k), m.postSamVec(1).nf(r), false, [], [], dispTrue, showToScreen, saveToFile);
	end;
elseif strcmp(par, 'lambda')
	resCell = cell(m.nr, 1+length(quantiles)+dispTrue);
	for r=1:m.nr
		if m.factorCov(r)
			for k=1:m.ncr(r)
				res = m.getPostLambda(r, k);
				label = strcat( 'lambda_', int2str(r), '_', int2str(k) );
				m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(r), m.ns, true, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
			end
		else
			res = m.getPostLambda(r, 0);
			label = strcat( 'lambda_', int2str(r) );
			resCell(r,:) = m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(r), m.ns, false, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
		end
	end
elseif strcmp(par, 'omega')
	resCell = cell(m.nr, 1+length(quantiles)+dispTrue);
	for r=1:m.nr
		if m.factorCov(r)
			warning('HMSC: cannot summarize omega for covariable-dependent level')
		else
			res = m.getPostOmega(r, []);
			label = strcat( 'omega_', int2str(r) );
			resCell(r,:) = m.summarize(label, res{1}, res{2}, quantiles, m.ns, m.ns, false, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
		end
	end
end

if ismember(par, {'beta', 'gamma', 'V', 'sigma', 'rho', 'ph', 'etas', 'lambdas', 'Am'} )
	resMean = resCell{1};
	varargout = resCell(2:end);
elseif ismember(par, {'alpha', 'eta', 'lambda', 'omega'} )
	resMean = resCell(:,1)';
	varargout = cell(1, nargout-1);
	for i=1:(size(resCell, 2)-1)
		varargout{i} = resCell(:,i)';
	end
end

end
