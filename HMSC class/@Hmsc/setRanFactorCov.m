function setRanFactorCov(m, XrCell)

m.factorCovScaleFlag = cell(1, m.nr);
m.factorCovScale = cell(1, m.nr);
if size(XrCell, 2) ~= m.nr
	error('HMSC: Size of Xr does not match number of random effects');
end
m.ncr = zeros(1, m.nr);
for r = 1:m.nr

% 	if( m.factorCov(r)==1 )
% 		if size(Xr{r}, 1) ~= m.np(r)
% 			error( strcat('HMSC: Dimension of levels for random effect', [' ' int2str(r)], ' does not match') )
% 		end
% 		m.ncr(r) = size(Xr{r}, 2);
% 	elseif m.factorCov(r)==2
% 		if size(Xr{r}, 2) ~= m.ns
% 			error( strcat('HMSC: Dimension species-specific random effect covariates at level', [' ' int2str(r)], ' must be specified as cell aray with length equal to number of species') )
% 		end
% 		
% 		m.ncr(r) = size(Xr{r}{1}, 2);
% 	else
% 		% possible to explicitely define Xr{r} as column of ones
% 		m.ncr(r) = 1;
% 	end
% 	m.factorCovScaleFlag{r} = zeros(1, m.ncr(r));
% 	m.factorCovScale{r} = zeros(2, m.ncr(r));
% 
% end
% m.Xr = Xr;

Xr = cell(1, m.nr);
ncr = nan(1, m.nr);
for r = 1:m.nr
	if( ~m.factorCov(r) && ~isempty(Xr{r}) )
		error( strcat('HMSC: Model is defined with non-covariating factors on effect', [' ' int2str(r)], ', but non-empty covariates matrix was passed for this effect') )
	end
	if m.factorCov(r)
		piKey = m.piKey{r}; 
		XrKey = XrCell{r}(:,1);
		XrMap = containers.Map(XrKey,1:length(XrKey));
		if ~all(XrMap.isKey(piKey))
			error('HMSC: some units, defined at level %d were not given factor covariates\n', r);
		end
		ind = cell2mat( XrMap.values(piKey) );
		Xr1 = XrCell{r}(ind, 2:size(XrCell{r}, 2));
		Xr1 = cell2mat(Xr1);
		Xr{r} = Xr1;
		ncr(r) = size(Xr{r}, 2);
	else
		ncr(r) = 1;
	end
	if( ~m.factorCov(r) && ~isempty(Xr{r}) )
		error( strcat('HMSC: Model is defined with non cavarying factors on effect', [' ' int2str(r)], ', but non-empty covariates matrix was passed for this effect') )
	end
end

m.XrCell = XrCell;
m.Xr = Xr;
m.ncr = ncr;
m.factorCovScaleFlag = cell(1, m.nr);
for r = 1:m.nr
	m.factorCovScaleFlag{r} = zeros(1, m.ncr(r));
end

end