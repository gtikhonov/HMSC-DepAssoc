function setDist(m, dist)

if strcmp(dist, 'normal')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 1; 
	m.dist(:, 2) = 1;
elseif strcmp(dist, 'probit')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 2; 
	m.dist(:, 2) = 0;
elseif strcmp(dist, 'poisson')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 3; 
	m.dist(:, 2) = 0;
elseif strcmp(dist, 'overdispersedPoisson')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 3; 
	m.dist(:, 2) = 1;
elseif strcmp(dist, 'tobitPoisson')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 4; 
	m.dist(:, 2) = 0;
elseif strcmp(dist, 'tobitOverdispersedPoisson')
	m.dist = zeros(m.ns, 4);
	m.dist(:, 1) = 4; 
	m.dist(:, 2) = 1;
else
	if size(dist,1) ~= m.ns
		error('HMSC: Distribution matrix number of rows must be equal to number of species');
	end
	if size(dist,2) ~= 2
		error('HMSC: Distribution matrix number of columns must be equal to 2');
	end
	dist = [dist, zeros(size(dist))];
	% todo - check correctness of dist
	m.dist = dist;
end

end