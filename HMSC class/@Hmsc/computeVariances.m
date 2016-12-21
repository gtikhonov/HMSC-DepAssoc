function [fixed, fixedsplit, random, traitR2] = computeVariances(m,group)

ngroups=max(group);
fixed=zeros(m.ns,1);
fixedsplit=zeros(m.ns,ngroups);
random=zeros(m.ns,m.nr);
traitR2=0;
cM=cov(m.X);

for i = 1:m.postSamN
	fixed1=zeros(m.ns,1);
	fixedsplit1=zeros(m.ns,ngroups);
	random1=zeros(m.ns,m.nr);
	beta=m.postSamVec(i).beta;
	if(m.traits)
		gamma=m.postSamVec(i).gamma;
		mu=((m.T)*gamma)';
		f=(m.X)*beta;
		a=(m.X)*mu;
		res=zeros(2,1);
		for j=1:m.ny
			tmp=cov(a(j,:),f(j,:));
			res(1)=res(1)+tmp(1,2)*tmp(2,1);
			res(2)=res(2)+tmp(1,1)*tmp(2,2);
		end
		traitR2=traitR2+res(1)/res(2);
	end
	for j=1:m.ns
		ftotal=beta(:,j)'*cM*beta(:,j);
		fixed1(j)=fixed1(j)+ftotal;
		for k=1:ngroups
			sel=(group==k);
			fpart=beta(sel,j)'*cM(sel,sel)*beta(sel,j);
			fixedsplit1(j,k)=fixedsplit1(j,k)+fpart;
		end
	end
	
	for level=1:m.nr
		lambda = m.postSamVec(i).lambda{level};
		nf=size(lambda,1);
		if m.factorCov(level)
			me=mean(m.Xr{level});
			cMr=cov(m.Xr{level});
			tM=cMr+me'*me;
			
			for factor =1:nf
				for k1=1:m.ncr(level)
					for k2=1:m.ncr(level)
						random1(:,level)=random1(:,level) + tM(k1,k2).*lambda(factor,:,k1)'.*lambda(factor,:,k2)';
					end
				end
			end
		else
			for factor =1:nf
				random1(:,level)=random1(:,level) + lambda(factor,:)'.*lambda(factor,:)';
			end
		end
	end
	if m.nr>0
		tot=sum([fixed1 random1]')';
		fixed=fixed+fixed1./tot;
		random=random+random1./repmat(tot,1,m.nr);
	else
		fixed=fixed+ones(m.ns,1);
	end
	fixedsplit=fixedsplit+fixedsplit1./repmat(sum(fixedsplit1')',1,ngroups);
end

fixed=fixed/m.postSamN;
random=random/m.postSamN;
fixedsplit=fixedsplit/m.postSamN;
traitR2=traitR2/m.postSamN;
