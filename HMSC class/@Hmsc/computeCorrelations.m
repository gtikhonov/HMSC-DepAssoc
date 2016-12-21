function [correlations,support,index] = computeCorrelations(m,level,Xr,threshold)

AOmega = m.getPostOmega(level, Xr);
npost=size(AOmega{1},1);
mR=zeros(m.ns,m.ns);
Rsig=zeros(m.ns,m.ns);
for i=1:npost
	Omega=reshape(AOmega{1}(i,:),m.ns,m.ns);
	[~,R]=cov2corr(Omega);
	mR=mR+R;
	Rsig=Rsig+(R>0);
end
Rsig=max(Rsig,npost-Rsig);
correlations=mR/npost;
support=Rsig/npost;
for i=1:m.ns
	support(i,i)=0;
end

index=1:m.ns;
scorrelations=correlations.*(support>threshold);
if(max(scorrelations(:))>0)
	[V,~]=eigs(scorrelations,2);
	angle=atan2(V(:,2),V(:,1));
	[~,index]=sort(angle);
	[~,ind] = min(V(index,2).^2 + V(index,1).^2);
	if ind > 1
		index = index([ind:m.ns, 1:(ind-1)]);
	end
end
