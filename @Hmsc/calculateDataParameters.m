function [iQg,detQg,iWgA,detWgA] = calculateDataParameters(m)

iQg = [];
detQg = [];
iWgA = [];
detWgA = [];
if m.phylogeny
	rhopw = m.rhopw;
	C = m.C;
	iC = corrcov(inv(m.C));
	for rg=1:size(rhopw, 1)
		rho = rhopw(rg, 1);
		if rho >= 0
			rhoC = rho*C;
		else
			rhoC = (-rho)*iC;
		end
		Q = rhoC+(1-abs(rho))*eye(m.ns);
		iQg = cat(3,iQg,inv(Q));
		cQ = chol(Q);
		detQg = [detQg; 2*sum(log(diag(cQ)))];
	end
end

for i=1:m.nr
	iWg=[];
	detWg=[];
	alphapw=m.alphapw{i};
	if m.spatial(i)
		xy = m.xy{i};
		di = zeros( m.np(i), m.np(i) );
		for j = 1:m.spatDim(i)
			xx = repmat(xy(:,j),1,m.np(i));
			dx = xx-xx';
			di = di+dx.^2;
		end
		distance = sqrt(di);
		for ag=1:size(alphapw, 1)
			alpha = alphapw(ag,1);
			if alpha < 1e-5
				W = eye(length(xy));
			else
				W = exp(-distance/alpha);
			end
			iWg = cat( 3, iWg, inv(W) );
			cholW = chol(W);
			detWg = [detWg; 2*sum(log(diag(cholW)))];
		end
	end
	iWgA{i} = iWg;
	detWgA{i} = detWg;
end
