function setCovScaling(m, scaleFlagVec)

if( length(scaleFlagVec)~= m.nc )
	error('HMSC: length of scaling flags vector must be equal to number of covariates');
end
if any(scaleFlagVec < 0) || any(scaleFlagVec > 2) || sum(scaleFlagVec==2)>1
	error('HMSC: elements of scaling flags could take integer values from 0 to 2 and only one could be intercept');
end
ind = scaleFlagVec==2;
if m.speciesX==false
	if any(ind) && any(m.X(:,ind)~=1)
		error('HMSC: intercept column must contain only ones');
	end
else
	for j=1:m.ns
		if any(ind) && any(m.X{j}(:,ind)~=1)
			error('HMSC: intercept column must contain only ones');
		end
	end
end

m.covScaleFlag = scaleFlagVec;
m.covScale = NaN(2, m.nc);

end

