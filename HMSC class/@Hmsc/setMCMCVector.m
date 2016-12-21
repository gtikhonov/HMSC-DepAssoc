function setMCMCVector(m, samVec)

if length(vector) < 1
	error('HMSC: Number of samples must be positive');
end
m.samples = NaN;
m.thinning = NaN;
m.samVec = samVec;

end