function setCovNames(m, covNames)

if length(covNames) ~= m.nc
	error('HMSC: Length of covariates names vector must be equal to number of covariates');
end;
m.covNames = covNames;

end