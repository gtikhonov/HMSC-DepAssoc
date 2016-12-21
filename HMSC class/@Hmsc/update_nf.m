function [nf,lambda,eta,alpha,psijh,delta] = update_nf(nf,ncr,ns,mcmc,np,nur,a2r,b2r,lambda,eta,alpha,psijh,delta,spatial)
nr = length(nf);
b0 = 1; b1 = 0.0005;
epsilon = 1e-3;                            % threshold limit
prop = 1.00;                               % proportion of redundant elements within columns
prob = 1/exp(b0 + b1*mcmc);                % probability of adapting
for i1=1:nr
	uu = rand;
	if uu < prob
		lambda1=lambda{i1};
		eta1=eta{i1};
		alpha1=alpha{i1};
		psijh1=psijh{i1};
		delta1=delta{i1};
		tmp=abs(lambda1)<epsilon;
		tmp = tmp(:,:)'; %stuck values for different covariates if lambda correspond to covariate-dependent factorA
		if ns>1
			tmp=sum(tmp);
		end
		lind = tmp/ns;    % proportion of elements in each column less than eps in magnitude
		vec = lind >=prop;
		num = sum(vec);       % number of redundant columns
		
		if mcmc > 20 && num == 0 && all(lind < 0.995)
			nf(i1) = nf(i1) + 1;
			lambda1(nf(i1),:,:) = 0;
			eta1(:,nf(i1))=normrnd(0,1,[np(i1),1]);
			alpha1(nf(i1)) = 1;
			psijh1(:,nf(i1),:) = gamrnd(nur/2,2/nur,[ns,1,ncr(i1)]);
			delta1(nf(i1),:) = gamrnd(a2r,1/b2r,[1,ncr(i1)]);
		elseif num > 0 && nf(i1)>2
			if num>nf(i1)-2
				nonred = [1 2];
			else
				nonred = setdiff(1:nf(i1),find(vec)); % non-redundant loadings columns
			end
			nf(i1) = length(nonred);
			lambda1 = lambda1(nonred,:,:);
			eta1 = eta1(:,nonred);
			if spatial(i1)
				alpha1 = alpha1(nonred);
			end
			psijh1 = psijh1(:,nonred,:);
			delta1 = delta1(nonred,:);
		end
		lambda{i1} = lambda1;
		eta{i1} = eta1;
		alpha{i1} = alpha1;
		psijh{i1} = psijh1;
		delta{i1} = delta1;
	end
end
