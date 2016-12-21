function setMCMCAdapt(m, adapt, fixNf, adaptXs, fixNfs)
if length(adapt) ~= 2
	error('HMSC: Adaptation parameter for latent factors must be a vector of 2 nonnegative integers');
end
%if ~islogical(fixNf)
%	error('HMSC: Fix number of factors parameter must be boolean');
%end
m.adapt = adapt;
m.fixNf = fixNf;

if m.includeXs
	if length(adaptXs) ~= 2
		error('HMSC: Adaptation parameter for dimesion reduction part must be a vector of 2 nonnegative integers');
	end
	if ~islogical(fixNfs)
		error('HMSC: Fix number of factors parameter must be boolean');
	end
	m.adaptXs = adaptXs;
	m.fixNfs = fixNfs;
end

end