function [nfs,lambdas,etas,psijhs,deltas] = update_nfs(mcmc,nus,a2s,b2s,lambdas,etas,psijhs,deltas)
[ny, nfs] = size(etas);
ns = size(lambdas, 2);

b0 = 1; b1 = 0.0005;
epsilon = 1e-3;                            % threshold limit
prop = 1.00;                               % proportion of redundant elements within columns
prob = 1/exp(b0 + b1*mcmc);                % probability of adapting

uu = rand;
if uu < prob
	lind = mean(abs(lambdas)<epsilon, 2);
	vec = lind>=prop;
	num = sum(vec);       % number of redundant columns
	
	if mcmc > 20 && num == 0 && all(lind < 0.995)
		nfs = nfs + 1;
		lambdas(nfs,:) = 0;
		etas(:,nfs)=normrnd(0,1,[ny,1]);
		psijhs(:,nfs) = gamrnd(nus/2,2/nus,[ns,1]);
		deltas(nfs) = gamrnd(a2s,1/b2s,1);
	elseif num > 0 && nfs>2
		if num > nfs-2
			nonred = [1 2];
		else
			nonred = setdiff(1:nfs, find(vec)); % non-redundant loadings columns
		end
		nfs = length(nonred);
		lambdas = lambdas(nonred,:);
		etas = etas(:,nonred);
		psijhs = psijhs(:,nonred);
		deltas = deltas(nonred);
	end
end

