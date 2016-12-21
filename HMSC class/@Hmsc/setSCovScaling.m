function setSCovScaling(m, scaleFlagVec)

if m.includeXs==0
	error('HMSC: model was defined without dimension reduced covariates Xs');
else
	if( length(scaleFlagVec)~= m.ncs )
		error('HMSC: length of dimension reduction scaling flags vector must be equal to number of dinesion reduced covariates');
	end
	if any(scaleFlagVec < 0) || any(scaleFlagVec > 2) || sum(scaleFlagVec==2)>1
		error('HMSC: elements of dimension reduction scaling flags could take integer values from 0 to 2 and only one could be intercept');
	end
	ind = scaleFlagVec==2;
	if m.includeXs==1
		if any(ind) && any(m.Xs(:,ind)~=1)
			error('HMSC: intercept column must contain only ones');
		end
	else
		for j=1:m.ns
			if any(ind) && any(m.Xs{j}(:,ind)~=1)
				error('HMSC: intercept column must contain only ones');
			end
		end
	end
	
	m.sCovScaleFlag = scaleFlagVec;
	m.sCovScale = NaN(2, m.nc);
end

end

