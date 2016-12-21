function setSpatialLocations(m, xyCell)

if size(xyCell, 2) ~= m.nr
	error('HMSC: Length of xy does not match number of levels of random effects');
end

xy = cell(1, m.nr);
spatDim = nan(1, m.nr);
for r = 1:m.nr
	if m.spatial(r)
		piKey = m.piKey{r}; 
		xyKey = xyCell{r}(:,1);
		xyMap = containers.Map(xyKey,1:length(xyKey));
		if ~all(xyMap.isKey(piKey))
			error('HMSC: some units, defined at spatially-defined level %d were not given spatial coordinates\n', r);
		end
		ind = cell2mat( xyMap.values(piKey) );
		xy1 = xyCell{r}(ind, 2:size(xyCell{r}, 2));
		xy1 = cell2mat(xy1);
		xy{r} = xy1;
		spatDim(r) = size(xy{r}, 2);
	end
	if( ~m.spatial(r) && ~isempty(xy{r}) )
		error( strcat('HMSC: Model is defined with non-spatial factors on level', [' ' int2str(r)], ', but non-empty locations matrix was passed for this effect') )
	end
end

m.xyCell = xyCell;
m.xy = xy;
m.spatDim = spatDim;

end