function setXs(m, Xs)

% check for correctness


if m.includeXs == 0 && ~isempty(Xs)
	error('HMSC: model was defined without dimension reduction part');
end


if m.includeXs == 2
	if length(Xs)~=m.ns
		error('HMSC: Length of Xs does not match mumber of species, defined by Y');
	end
	m.ncs = size(Xs{1}, 2);
	for i = 1:m.ns
		if any(size(Xs{i}) ~= [m.ny, m.ncs])
			error(paste('HMSC: Xs{', int2str(i), '} dimensions are not compatible with specified'));
		end
	end
else
	m.ncs = size(Xs, 2);
	if any(size(Xs) ~= [m.ny, m.ncs])
		error('HMSC: Xs dimensions are not compatible with specified');
	end
end
m.Xs = Xs;
m.sCovScaleFlag = zeros(1, m.ncs);
m.sCovScale = zeros(2, m.ncs);

end