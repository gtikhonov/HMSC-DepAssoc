function setFactorCovScaling(m, scaleFlagVecList)

for r=1:m.nr
	if m.factorCov(r) > 0
		scaleFlagVec = scaleFlagVecList{r};
		if( length(scaleFlagVec)~= m.ncr(r) )
			error('HMSC: length of factor scaling flags vector must be equal to number of factor covariates');
		end
		if any(scaleFlagVec < 0) || any(scaleFlagVec > 2) || sum(scaleFlagVec==2)>1
			error('HMSC: elements of dimension reduction scaling flags could take integer values from 0 to 2 and only one could be intercept');
		end
		ind = scaleFlagVec==2;
		if any(ind)
			if m.factorCov(r)==1 
				if any(m.Xr{r}(:,ind)~=1)
					error('HMSC: intercept column must contain only ones');
				end
			else
				for j=1:m.ns
					if any(m.Xr{r}{j}(:,ind)~=1)
						error('HMSC: intercept column must contain only ones');
					end
				end
			end
		end
		m.factorCovScaleFlag{r} = scaleFlagVec;
		m.factorCovScale{r} = NaN(2, m.nc);
	end
end


end

