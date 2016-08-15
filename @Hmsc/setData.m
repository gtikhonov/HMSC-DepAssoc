function setData(m, Y, dist, X, piCell, xyCell, T, C, XrCell, Xs)

[ny, ns] = size(Y);
if m.speciesX
	nc = size(X{1}, 2);
else
	nc = size(X, 2);
end
setDim(m, ny, ns, nc);
m.setY(Y);
m.setX(X);
m.setDist(dist);
m.setPi(piCell);
if ~isempty(T)
	setTraits(m, T);
end
if ~isempty(C)
	setPhylogeny(m, C);
end
if ~isempty(xyCell)
	setSpatialLocations(m, xyCell)
end
if ~isempty(XrCell)
	m.setRanFactorCov(XrCell)
end
if ~isempty(Xs)
	m.setXs(Xs)
end

end