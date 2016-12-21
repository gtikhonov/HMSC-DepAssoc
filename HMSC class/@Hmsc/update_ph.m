function ph = update_ph(T,beta,iV,ph,gamma,nu)
nc = length(gamma);
[ns nt]=size(T);
aph = (nu+nc)/2;
mbeta = gamma'*T';
for j = 1:ns
	bph = 0.5*(nu + (beta(:,j)-mbeta(:,j))'*iV*(beta(:,j)-mbeta(:,j)));
	ph(j) = gamrnd(aph,1/bph);
end
