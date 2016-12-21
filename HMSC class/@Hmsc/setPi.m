function setPi(m, piCell)

if size(piCell, 2) ~= m.nr
	error('HMSC: Number of columns of pi must be equal to number of levels')
end

np = nan(1, m.nr);
for r=1:m.nr
	np(r) = length(unique(piCell(:,r)));
end

m.piCell = piCell;
m.np = np;

pi = nan(m.ny, m.nr);
unCell = cell(1, m.nr);
mapCell = cell(1, m.nr);
for r=1:m.nr
	keySet = unique(piCell(:,r));
	if length(keySet)==m.ny
		keySet = piCell(:,r);
	end
	unCell{r} = keySet;
	valueSet = 1:length(keySet);
	mapCell{r} = containers.Map(keySet,valueSet);
	pi(:,r) = cell2mat(mapCell{r}.values( piCell(:,r)) );
end
m.pi = pi;
m.piKey = unCell;
m.piMap = mapCell;

end