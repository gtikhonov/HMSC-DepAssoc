function genXs(m, ncs)

if m.includeXs == 0
	error('HMSC: model was defined without dimension reduction part');
end

if m.includeXs == 2
	if isempty( m.ns )
		error('HMSC: Number of species must be predefined to generate species-specific Xs');
	end
	Xs = cell(1, m.ns);
	for i = 1:m.ns
		Xs1 = normrnd(0,1,[m.ny,ncs]);
		Xs{i}=Xs1;
	end
else
	Xs = normrnd(0,1,[m.ny,ncs]);
end
m.setXs(Xs);

sCovNames = cell(1,ncs);
for i=1:ncs
	sCovNames{i} = strcat('sCov_',int2str(i));
end
m.setSCovNames(sCovNames);

end