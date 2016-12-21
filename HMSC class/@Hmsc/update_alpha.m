function alpha = update_alpha(eta,alpha,nf,pi,spatial,iWg,iWgR,detWg,alphapw)
% eta - cell array of latent factors at different levels. Each cell corresponds to one level. Each cell contains one
% matrix with dimention of number_of_units_for_this_level to number_of_factors_on this level
% alpha - cell array of parameters of spatial autocorrelation on different levels. Each cell corresponds to one level.
% Within each cell is either an empty object for non-spatial levels, or an array of indices for spatial level. These
% indices are combined with discrite grid for alpha for that level, given by alphapw{level}(:,1)
% nf - array with number of factors at different levels
% pi - matrix of correspondence of observations to units at different levels
% spatial - binary array, indicating whether the level is spatial or not
% iWg - cell array of precalculated inverse matrices
% detWg - cell array of precalculated determinants
% alphapw - cell array of discrite grids for alpha at different levels and the prior probabilities of those 
if min(pi)<0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
end
% iterate over the latent factors levels, which are defined as spatial
for r=1:nr
	%user specified: layers
	if spatial(r)==true
% If we look at one particular spatial level, Wg is the exp(-alpha*DIST_MATRIX), its inverse and discriminant are 
% precalculated for every alpha from discrete grid (locations give by alphapw{r}(:,1), prior weights by alphapw{r}(:,1))
% just prior to the beginning of sampling.
		iWg1 = iWg{r};
      iWgR1 = iWgR{r};
		detWg1 = detWg{r};
		alphapw1 = alphapw{r};
		alpha1=alpha{r};
		eta1=eta{r};
		gridN = size(alphapw1, 1);
% iterate over the random factors at this spatial level (for each alpha is sampled independently on any other)
      tmpMat = permute(sum(mtimesx(iWgR1, eta1).^2), [3,2,1]);
%       tmpMat = nan(gridN, nf(r));
%       for ag = 1:gridN
%          a = iWgR1(:,:,ag)*eta1;
%          tmpMat(ag,:) = sum(a.^2);
%       end
		for h=1:nf(r);
% 			tmp = nan(gridN, 1);
% 			etaVec = eta1(:,h);
% Next 4 lines calculates the posterior log-likelihood for different discrete values of alpha - the complete description of 
% why this formula is correct is given in the doc file with derivations of sampling formulas.
% Basicly it comes from the fact the etaVec ~ N([0,...0], Wg), so we just calculate the likelihoods for different Wg
% (which are in turn defined by different alpha from the discrete grid)
%          etaVecT = etaVec';
%          for ag = 1:gridN
% 				tmp(ag) = etaVecT*iWg1(:,:,ag)*etaVec;
%          end
%          for ag = 1:gridN
%             a = iWgR1(:,:,ag)*etaVec;
%             tmp(ag) = sum(a.^2);
%          end
%          tmp = sum(sum(bsxfun(@times,permute(iWgR1, [1,3,2]),permute(etaVec, [2,3,1])), 3).^2)';
         tmp = tmpMat(:,h);
			like = log(alphapw1(:,2))-1/2*detWg1-(1/2)*tmp;
% Recenter the log-likelihood to be from -inf to 0
			like=like-max(like);
% Convert to likelihood
			like = exp(like);
% Normalize the likelihood
			like = like/sum(like);
% sample alpha for given random factor on given spatial level from resulted likelihood
			alpha1(h)=randsample(length(like),1,true,like);
		end
		alpha{r}=alpha1;
	end
end

