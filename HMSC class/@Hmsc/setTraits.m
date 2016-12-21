function setTraits(m, T)

if ~m.traits
	error('HMSC: Model was defined without traits');
end
nsIn = size(T, 1);
if nsIn ~= m.ns
	error('HMSC: Number of rows of T should be equal to number of species');
end
m.T = T;
m.nt = size(m.T, 2);

end