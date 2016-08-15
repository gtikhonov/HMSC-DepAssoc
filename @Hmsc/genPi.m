function genPi(m, np)
if length(np) ~= m.nr
	error('HMSC: Length of number of levels for random factors must be equal to the number of latent factors');
end

if m.nr>0
	pi=nan(m.ny, m.nr);
	for r = 1:m.nr
		pi(:,r) = rem((1:m.ny)-1, np(r))+1;
		if r > 1
			pi(:,r) = pi(randperm(m.ny), r);
		end
	end
end
pi = num2cell(pi);
pi = cellfun(@num2str, pi, 'UniformOutput', false);
for r=1:m.nr
	pi(:,r) = cellfun(@(c) [num2str(r), '-', c], pi(:,r), 'UniformOutput', false);
end
m.setPi(pi)

end