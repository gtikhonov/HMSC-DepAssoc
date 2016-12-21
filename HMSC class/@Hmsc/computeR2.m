function R2 = computeR2(m,predN)

predList = m.predict(predN, m.X, m.piCell, m.xyCell, m.XrCell, true);
pred = predList{1};
for i=2:predN
	pred = pred+predList{i};
end
pred = pred / predN;

samplingUnitIncluded=false;
for level=1:m.nr
	if max(m.pi(:,level))==m.ny
		samplingUnitIncluded=true;
	end
end
z=m.nr+(~samplingUnitIncluded);
R2 = repmat(NaN,m.ns,z);
for level=1:z
	if level<=m.nr
		nunits=max(m.pi(:,level));
	else
		nunits=0;
	end
	if or(nunits==m.ny,level>m.nr)
		for j=1:m.ns
			if m.dist(j)==2 %probit model
				R2(j,level) = mean(pred(m.Y(:,j)==1,j))-mean(pred(m.Y(:,j)==0,j));
			else  %Poisson or normal models
				sel=~isnan(m.Y(:,j));
				if sum(sel)>4
					R2(j,level) = corr(pred(sel,j),m.Y(sel,j));
				end
			end
		end
	else
		obs=repmat(NaN,nunits,m.ns);
		preds=repmat(NaN,nunits,m.ns);
		for j=1:m.ns
			for i=1:nunits
				sel=and(~isnan(m.Y(:,j)),m.pi(:,level)==i);
				if(sum(sel)>0)
					obs(i,j)=sum(m.Y(sel,j));
					preds(i,j)=sum(pred(sel,j));
				end
			end
		end
		for j=1:m.ns
			sel=~isnan(obs(:,j));
			if sum(sel)>4
				R2(j,level) = corr(obs(sel,j),preds(sel,j));
			end
		end
	end
end

