function setMCMCOptions(m, samples, thinning)
if samples < 1
	error('HMSC: Number of samples must be positive');
end
if thinning < 1
	error('HMSC: Thinning value must be positive');
end

m.samples = samples;
m.thinning = thinning;
m.samVec = (1:samples) * thinning;

end