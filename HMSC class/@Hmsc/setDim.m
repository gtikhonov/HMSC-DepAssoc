function setDim(m, ny, ns, nc)

m.ny = ny;
m.ns = ns;
m.nc = nc;
if ~m.traits
	m.T = ones(m.ns, 1);
	m.nt = 1;
end

end