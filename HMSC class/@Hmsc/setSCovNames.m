function setSCovNames(m, sCovNames)

if length(sCovNames) ~= m.ncs
	error('HMSC: Length of dimension-reduced covariates names vector must be equal to number of dimension-reduced covariates');
end;
m.sCovNames = sCovNames;

end