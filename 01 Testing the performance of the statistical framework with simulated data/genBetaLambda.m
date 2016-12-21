function genBetaLambda(globalCount, commType, assocType, nf, folder)
rng(globalCount);
nc = 2;
ns = 50;
ncr = 2;
betaT = nan(nc, ns);
if strcmp(commType, 'poor')
   betaT(1,:) = -2+randn(1, ns);
elseif  strcmp(commType,'rich')
   betaT(1,:) = -1+randn(1, ns);
else
   error('Wrong community type')
end
betaT(2:nc,:) = randn(nc-1, ns);

lambdaT = nan(nf,ns,ncr);
mult = 0.4/sqrt(nf/2);
aaa
lambdaT(:,:,1) = randn(nf, ns);
switch assocType
   case 'const'
      lambdaT(:,:,2) = 0;
   case 'center'
      lambdaT(:,:,1) = lambdaT(:,:,1) + (0.9 + 0.1*randn(nf, ns, ncr-1));
      lambdaT(:,:,2) = 0;
   case 'stress'
      lambdaT(:,:,2:ncr) = 0.9 + 0.1*randn(nf, ns, ncr-1);
      lambdaT(:,:,1) = lambdaT(:,:,1) + lambdaT(:,:,2);
   case 'random'
      lambdaT(:,:,2:ncr) = 1.8*randn(nf, ns, ncr-1);
   case 'equal'
      lambdaT = 1.25*randn(nf, ns, ncr);
   otherwise
      error('Wrong association type') 
end
lambdaT = mult*lambdaT;


save(fullfile(folder, sprintf('betaLambda %.3d.mat', globalCount)), 'betaT', 'lambdaT')
