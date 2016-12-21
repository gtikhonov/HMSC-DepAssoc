function genTraits(m, nt)
if nt < 1
	error('HMSC: Number of traits must be a positive value');
end

T = normrnd(0, 1, [m.ns, nt]);
T(:, 1) = 1;
m.setTraits(T);

traitNames = cell(1, nt);
for i=1:nt
	traitNames{i} = strcat('trait_', int2str(i));
end
m.setTraitNames(traitNames);

end