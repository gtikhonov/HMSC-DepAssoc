function genSpatialLocations(m, spatDim)
if length(spatDim) ~= m.nr
	error('HMSC: Length of vector with spatial dimentions for random facors must be eaual to number of random factors');
end

xy = cell(1, m.nr);

for r=1:m.nr
	if m.spatial(r)
		xy1=rand(m.np(r), spatDim(r));
		xy{r} = [m.pi(:,r), num2cell(xy1)];
	end
end

m.setSpatialLocations(xyCell)

end