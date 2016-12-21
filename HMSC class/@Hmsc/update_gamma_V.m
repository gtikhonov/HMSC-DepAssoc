function [gamma,iV] = update_gamma_V(T,beta,gamma,ph,rho,V0,f0,Ugamma,mgamma,phylogeny,iQg,outlierspecies)
[ns nt]=size(T);
[nc ns]=size(beta);

if phylogeny
	iQ=iQg(:,:,rho);
else
	iQ=eye(ns);
end

Ts = T;
betas=beta;
if outlierspecies
	for i=1:ns
		Ts(i,:)=sqrt(ph(i))*T(i,:);
		betas(:,i)=sqrt(ph(i))*beta(:,i);
	end
end

E=betas'-(Ts*gamma);
A=E'*iQ*E;
Vn=inv(A+V0);
Vn=(Vn+Vn')/2;
iV = wishrnd(Vn,f0+ns);

tmp=Ts'*iQ;
Ugammas=inv(inv(Ugamma)+kron(iV,tmp*Ts));
Ugammas=(Ugammas+Ugammas')/2;
res=tmp*betas'*iV;
mgammas=Ugammas*(inv(Ugamma)*mgamma+res(:));
gamma=mvnrnd(mgammas,Ugammas)';
gamma=reshape(gamma,nt,nc);
