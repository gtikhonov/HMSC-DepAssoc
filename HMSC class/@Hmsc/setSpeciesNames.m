function setSpeciesNames(m, spNames)

if length(spNames) ~= m.ns
	error('HMSC: Length of species names vector must be equal to number of species');
end;
m.spNames = spNames;

end