function createPostSamVec(m, runs, thinVec)
ind = [];
vec = [];

for i = 1:length(runs)
	k = runs(i);
	shift = (k-1) * m.samVec(end);
	flag = 0;
	if isempty(m.repPar{k})
		m.postFileToRam(k);
		flag = 1;
	end;
	ind = [ind, shift + m.samVec(thinVec{i})];
	vec = [vec, m.repPar{k}{thinVec{i}}];
	if flag
		m.repPar{k} = [];
	end;
end

m.postSamInd = ind;
m.postSamVec = vec;
m.postSamN = length(ind);

end