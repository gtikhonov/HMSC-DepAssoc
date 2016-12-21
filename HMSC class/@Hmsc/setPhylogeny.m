function setPhylogeny(m, C)

if ~m.phylogeny
	error('HMSC: Model was defined without phylogeny');
end
if( size(C, 1) ~= m.ns || size(C, 2) ~= m.ns )
	error('HMSC: C must be a square matrix with nunber of rows and columns equal to number of species');
end
if C ~= transpose(C)
	error('HMSC: C must be a symmetric matrix');
else
	% check if positively definite
	[~,p] = chol(C);
	if p > 0
		error('HMSC: C must be a positive definite matrix');
	end
end
m.C = C;

end