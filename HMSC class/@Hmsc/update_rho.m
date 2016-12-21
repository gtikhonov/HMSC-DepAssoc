function rho = update_rho(rho,beta,T,gamma,detQg,iQg,iV,rhopw,ph,outlierspecies)
% rho - index of rho on discrete grid (rhopw)
% beta - responses of species to covariates
% T - matrix of species traits
% gamma - matrix of impacts of Traints on beta
% detQg - precalculated determinants of Qg for discrete grid of rho
% iQg - precalculated inverses of Qg for discrete grid of rho
% iV - inverse of parameter matrix in inv-Wishart distribution, used as prior for covariation matrix of beta atop of the part, explained by traits
% rhopw - discrete grid for rho and the prior probabilities of those points
% ph - array of multipliers, used for outlier species overdispersion
% outlierspecies - flag, whether outlier species are present in this model
[nc ns]=size(beta);
tmp=[];
% subtract the part of responses, which is explained by traits. Qg = rho*C + (1-rho)*I (it is like this only for
% positive rho, but should not matter now). Here C is the user-defined phylogenic correlation matrix, Qg - corresponding
% corelation matrix for responses. Then basically we 
res=beta'-(T*gamma);
res=res(:);
% detQg - precalculated determinants of
% thus this cycle just iterates over the discrete grid of rho, which locations are given by rhopw(:,1)
% and prior weights for those locations are given by rhopw(:,2) 
for rg = 1:length(detQg)
	if outlierspecies
% make the correction for outliers
		PHI=(1./sqrt(ph))*(1./sqrt(ph))';
		iVs = kron(iV,(1./PHI).*iQg(:,:,rg));
	else
% This holds as each column of beta'-(T*gamma) is supposed to be independently distributed as N([0,...,0], Qg) thus here
% is the inverse of the covariance matrix for all columns stacked
		iVs = kron(iV,iQg(:,:,rg));
	end
% Just store the quadratic part of likelihood expression for different Qg (different rho)
	tmp = [tmp; res'*iVs*res];
end
% In next two lines the log-likelihood for different rho from the discrete grid is calculated. Basically it is just the
% multivariate normal distribution likelihood formula, adjusted by prior weights for rho, given by rhopw(:,2).
logdetg=-ns*log(det(iV))+nc*detQg;
like = log(rhopw(:,2))-1/2*logdetg-(1/2)*tmp;
% Recenter the log-likelihood to be from -inf to 0
like = like-max(like);
% Convert to likelihood
like = exp(like);
% Normalize the likelihood
like = like/sum(like);
rho=randsample(length(like),1,true,like);
