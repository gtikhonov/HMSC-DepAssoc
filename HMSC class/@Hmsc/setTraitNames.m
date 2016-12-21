function setTraitNames(m, traitNames)

if length(traitNames) ~= m.nt
	error('HMSC: Length of traits names vector must be equal to number of traits');
end;
m.traitNames = traitNames;

end