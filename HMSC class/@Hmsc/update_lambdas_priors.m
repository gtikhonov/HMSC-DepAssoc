function [psijhs,deltas] = update_lambdas_priors(nfs,nus,a1s,a2s,b1s,b2s,psijhs,deltas,lambdas)
ns = size(lambdas,2);
tauhs = cumprod(deltas);
psijhs = gamrnd(nus/2 + 0.5,1./(nus/2 + 0.5*bsxfun(@times,(lambdas').^2,tauhs')));
mat = bsxfun(@times,psijhs,(lambdas').^2);
ad = a1s + 0.5*ns*nfs;
bd = b1s + 0.5*(1/deltas(1))*sum(tauhs.*sum (mat)');
deltas(1) = gamrnd(ad,1/bd);
tauhs = cumprod(deltas);
for h = 2:nfs
	ad = a2s + 0.5*ns*(nfs-h+1);
	bd = b2s + 0.5*(1/deltas(h))*sum(tauhs (h:end).*sum (mat(:,h:end))');
	deltas(h) = gamrnd(ad,1/bd);
	tauhs = cumprod(deltas);
end
