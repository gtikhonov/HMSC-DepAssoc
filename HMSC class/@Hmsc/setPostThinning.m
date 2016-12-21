function setPostThinning(m, runs, postThin)

if isnan(m.samples)
	error('HMSC: MCMC sampling was not specified with equal thinning via setMCMCOptions method. Use createPostSamVec method to sample posterior with given vector');
end
thinVec = find( rem( 1:m.samples, postThin ) == 0 );
thinVec = repmat({thinVec}, 1, length(runs) );

m.createPostSamVec(runs, thinVec)

end