function setX(m, X)

if m.speciesX
	if length(X)~=m.ns
		error('HMSC: Length of X does not match mumber of species, defined by Y');
	end
	for i = 1:m.ns
		if any(size(X{i}) ~= [m.ny, m.nc])
			error(paste('HMSC: X{', int2str(i), '} dimentions are not compatible with specified'));
		end
	end
else
	if any(size(X) ~= [m.ny, m.nc])
		error('HMSC: X dimentions are not compatible with specified');
	end
end
m.X = X;
m.covScaleFlag = zeros(1, m.nc);
m.covScale = zeros(2, m.nc);

end