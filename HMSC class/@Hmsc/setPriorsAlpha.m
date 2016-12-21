function setPriorsAlpha(m, alphap, alphaw)

if length(alphap) ~= m.nr
	error('HMSC: First parameter must be a cell array with length equal to number of levels of latent factors');
end
if length(alphaw) ~= m.nr
	error('HMSC: Second parameter must be a cell array with length equal to number of levels of latent factors');
end

m.alphapw = cell(1, m.nr);
for r=1:m.nr
	if m.spatial(r)
		if isempty(alphap{r}) && ~isempty(alphaw{r})
			error('HMSC: Descrete grid at level %d must be explicitely specified when explicitely setting the weights at level %d', r, r);
		end
		if ~isempty(alphap{r})
			ap = alphap{r};
			if any(ap<0)
				error('HMSC: Descrete grid for factor level %d must be positive', r);
			end
		else
			d = sqrt(sum((max(m.xy{r})-min(m.xy{r})).^2));
			alphap = (0:0.01:1)*d;
		end
		if ~isempty(alphaw{r})
			aw = alphaw{r};
			if length(aw) ~= length(ap)
				error('HMSC: Length of weights vector must be equal to the length of grid vector at level %d', r);
			end
			if any(aw) < 0
				error('HMSC: Weights for level %d must be positive', r);
			end
			if abs(sum(aw)-1) > 0.001
				error('HMSC: Sum of weights at level %d must be equal to 1', r);
			end
		else
			aw = ones(1,length(ap));
			aw = 0.5*aw/(sum(aw)-1);
			aw(1) = 0.5;
		end
		m.alphapw{r} = [ap; aw]';
	else
		if ~isempty(alphap{r})
			error('HMSC: Random factor level %d was specified as non-spatial, but a non empty grid was provided', r );
		end
		if ~isempty(alphaw{r})
			error('HMSC: Random factor level %d was specified as non-spatial, but non empty weights were provided', r );
		end
	end
end

end