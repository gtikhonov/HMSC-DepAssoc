function genRanFactorCov(m, ncr)
if length(ncr) ~=  m.nr
	error('HMSC: Length of vector with number of covariates for latent factors must be equal to number of latent factors');
end
ncr( isnan(ncr) ) = 1;

Xr = cell(1, m.nr);
for r = 1:m.nr
	if m.factorCov(r) == 1
		Xr1=normrnd(0, 1, [m.np(r), ncr(r)]);
		Xr1(:, 1) = 1;
		Xr{r} = [m.pi(r,:), num2cell(Xr1)];
	elseif m.factorCov(r) == 2
		Xr{r} = cell(1, m.ns);
		for j = 1:m.ns
			Xr1=normrnd(0, 1, [m.np(r), ncr(r)]);
			Xr1(:, 1) = 1;
			Xr{r}{j} = [m.pi(r,:), num2cell(Xr1)];
		end
	else
		Xr{r} = [];
	end
end

m.setRanFactorCov(XrCell);
end