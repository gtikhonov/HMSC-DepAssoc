function summary(m, par, quantiles, dispTrue, showToScreen, saveToFile)

% check if can dispTrue

if strcmp(par, 'beta')
	res = m.getPostBeta();
	m.summarize('beta', res{1}, res{2}, quantiles, m.nc, m.ns, true, m.spNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'gamma')
	res = m.getPostGamma();
	m.summarize('gamma', res{1}, res{2}, quantiles, m.nt, m.nc, false, m.traitNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'V')
	res = m.getPostV();
	m.summarize('V', res{1}, res{2}, quantiles, m.nc, m.nc, false, m.covNames, m.covNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'sigma')
	res = m.getPostSigma();
	m.summarize('sigma', res{1}, res{2}, quantiles, m.ns, 1, false, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'rho')
	res = m.getPostRho();
	m.summarize('rho', res{1}, res{2}, quantiles, 1, 1, false, {'rho'}, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'ph')
	res = m.getPostPh();
	m.summarize('ph', res{1}, res{2}, quantiles, m.ns, 1, false, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'etas')
	res = m.getPostEtas();
	m.summarize('etas', res{1}, res{2}, quantiles, m.ncs, m.postSamVec(1).nfs, false, m.sCovNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'lambdas')
	res = m.getPostLambdas();
	m.summarize('lambdas', res{1}, res{2}, quantiles, m.postSamVec(1).nfs, m.ns, true, m.spNames, [], dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'Am')
	res = m.getPostAm();
	m.summarize('Am', res{1}, res{2}, quantiles, m.ncs, m.ns, true, m.spNames, m.sCovNames, dispTrue, showToScreen, saveToFile);
elseif strcmp(par, 'alpha')
	for i=1:m.nr
		if m.spatial(i)
			res = m.getPostAlpha(i);
			label = strcat( 'alpha_', int2str(i) );
			m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(i), 1, false, [], [], dispTrue, showToScreen, saveToFile);
		end
	end
elseif strcmp(par, 'eta')
	for k = 1:m.nr
		res = m.getPostEta(k);
		label = strcat('eta_', int2str(k));
		m.summarize(label, res{1}, res{2}, quantiles, m.np(k), m.postSamVec(1).nf(i), false, [], [], dispTrue, showToScreen, saveToFile);
	end;
elseif strcmp(par, 'lambda')
	for i=1:m.nr
		if m.factorCov(i)
			for k=1:m.ncr(i)
				res = m.getPostLambda(i, k);
				label = strcat( 'lambda_', int2str(i), '_', int2str(k) );
				m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(i), m.ns, true, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
			end
		else
			res = m.getPostLambda(i, 0);
			label = strcat( 'lambda_', int2str(i) );
			m.summarize(label, res{1}, res{2}, quantiles, m.postSamVec(1).nf(i), m.ns, false, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
		end
	end
elseif strcmp(par, 'omega')
	for i=1:m.nr
		if m.factorCov(i)
		else
			res = m.getPostOmega(i, []);
			label = strcat( 'omega_', int2str(i) );
			m.summarize(label, res{1}, res{2}, quantiles, m.ns, m.ns, false, m.spNames, m.spNames, dispTrue, showToScreen, saveToFile);
		end
	end
end

end
