function genX(m)
if isempty( m.nc )
	error('HMSC: Number of covariates must be predefined to generate X');
end
if isempty( m.ny )
	error('HMSC: Number of sampling units must be predefined to generate X');
end
if m.speciesX
	if isempty( m.ns )
		error('HMSC: Number of species must be predefined to generate species-specific X');
	end
	X = cell(1, m.ns);
	for i = 1:m.ns
		X1 = normrnd(0,1,[m.ny,m.nc]);
		X1(:,1) = 1;
		X{i}=X1;
	end
else
	X = normrnd(0,1,[m.ny,m.nc]);       % habitat covariates
	X(:,1) = 1;
end
m.setX(X);

covNames = cell(1,m.nc);
for i=1:m.nc
	covNames{i} = strcat('cov_',int2str(i));
end
m.setCovNames(covNames);

spNames = cell(1,m.ns);
for i=1:m.ns
	spNames{i} = strcat('sp_',int2str(i));
end
m.setSpeciesNames(spNames);
end