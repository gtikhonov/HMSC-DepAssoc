classdef HmscPar
	properties
		sigma, gamma, beta, V
		nf, lambda, eta, delta, psijh
		alpha % spatial random factors
		etas, lambdas, nfs, deltas, psijhs % includeXs
		rho % phylogeny
		ph % outlier species
	end
	methods
		function writeToText(p, path)
			fileID = fopen(path, 'w');
			
			fprintf(fileID, 'sigma ');
			if isempty(p.sigma)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				ns = size(p.sigma, 1);
				fprintf(fileID,'%d\n', ns);
				fprintf(fileID,'%f ', diag(p.sigma));
				fprintf(fileID,'\n');
			end
			
			fprintf(fileID, 'gamma ');
			if isempty(p.gamma)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				[nt, nc] = size(p.gamma);
				fprintf(fileID,'%d %d\n', nt, nc);
				fprintf(fileID,'%f ', p.gamma);
				fprintf(fileID,'\n');
			end
			
			fprintf(fileID, 'beta ');
			if isempty(p.beta)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				[nc, ns] = size(p.beta);
				fprintf(fileID,'%d %d\n', nc, ns);
				fprintf(fileID,'%f ', p.beta);
				fprintf(fileID,'\n');
			end
			
			fprintf(fileID, 'V ');
			if isempty(p.V)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nc = size(p.V, 1);
				fprintf(fileID,'%d\n', nc);
				fprintf(fileID,'%f ', p.V);
				fprintf(fileID,'\n');
			end
			
			fprintf(fileID, 'nf ');
			if isempty(p.nf)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.nf);
				fprintf(fileID,'%d\n', nr);
				fprintf(fileID,'%d ', p.nf);
				fprintf(fileID,'\n');
			end
			
			fprintf(fileID, 'lambda ');
			if isempty(p.lambda)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.lambda);
				fprintf(fileID,'%d\n', nr);
				for r = 1:nr
					[nf, ns] = size(p.lambda{r});
					fprintf(fileID,'%d %d\n', nf, ns);
					fprintf(fileID,'%f ', p.lambda{r});
					fprintf(fileID,'\n');
				end
			end
			
			fprintf(fileID, 'eta ');
			if isempty(p.eta)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.eta);
				fprintf(fileID,'%d\n', nr);
				for r = 1:nr
					[np, nf] = size(p.eta{r});
					fprintf(fileID,'%d %d\n', np, nf);
					fprintf(fileID,'%f ', p.eta{r});
					fprintf(fileID,'\n');
				end
			end
			
			fprintf(fileID, 'delta ');
			if isempty(p.delta)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.delta);
				fprintf(fileID,'%d\n', nr);
				for r = 1:nr
					nf = length(p.delta{r});
					fprintf(fileID,'%d\n', nf);
					fprintf(fileID,'%f ', p.delta{r});
					fprintf(fileID,'\n');
				end
			end
			
			fprintf(fileID, 'psijh ');
			if isempty(p.psijh)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.psijh);
				fprintf(fileID,'%d\n', nr);
				for r = 1:nr
					[ns, nf] = size(p.psijh{r});
					fprintf(fileID,'%d %d\n', ns, nf);
					fprintf(fileID,'%f ', p.psijh{r});
					fprintf(fileID,'\n');
				end
			end
			
			fprintf(fileID, 'alpha ');
			if isempty(p.alpha)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				nr = length(p.alpha);
				fprintf(fileID,'%d\n', nr);
				for r = 1:nr
					if ~isnan(p.alpha{r})
						nf = length(p.alpha{r});
						fprintf(fileID,'%d\n', nf);
						fprintf(fileID,'%d ', p.alpha{r});
					else
						fprintf(fileID,'-1');
					end
					fprintf(fileID,'\n');
				end
			end
			
			fprintf(fileID, 'rho ');
			if isempty(p.rho)
				fprintf(fileID,'0\n');
			else
				fprintf(fileID,'1\n');
				if ~isnan(p.rho)
					fprintf(fileID,'%d', p.rho);
				else
					fprintf(fileID,'-1');
				end
				fprintf(fileID,'\n');
			end
			
			fclose(fileID);
		end
	end
	methods (Static = true)
		function p = readFromText(path)
			fileID = fopen(path, 'r');
			sigma = [];
			gamma = [];
			beta = [];
			V = [];
			nf = [];
			lambda = [];
			eta = [];
			delta = [];
			psijh = [];
			alpha = [];
			rho = [];
			
			ex = fscanf(fileID, 'sigma %d');
			if ex
				ns = fscanf(fileID, '%d', 1);
				sigma = fscanf(fileID, '%f', ns)';
			end
			
			ex = fscanf(fileID, '\ngamma %d');
			if ex
				s = fscanf(fileID, '%d %d', 2)';
				gamma = fscanf(fileID, '%f ', s);
			end
			
			ex = fscanf(fileID, '\nbeta %d');
			if ex
				s = fscanf(fileID, '%d %d', 2)';
				beta = fscanf(fileID, '%f ', s);
			end
			
			ex = fscanf(fileID, '\nV %d');
			if ex
				nt = fscanf(fileID, '%d', 1);
				V = fscanf(fileID, '%f ', [nt, nt]);
			end
			
			ex = fscanf(fileID, '\nnf %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				nf = fscanf(fileID, '%d ', nr)';
			end
			
			ex = fscanf(fileID, '\nlambda %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				lambda = cell(1, nr);
				for r = 1:nr
					s = fscanf(fileID, '%d %d', 2)';
					lambda{r} = fscanf(fileID, '%f ', s);
				end
			end
			
			ex = fscanf(fileID, '\neta %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				eta = cell(1, nr);
				for r = 1:nr
					s = fscanf(fileID, '%d %d', 2)';
					eta{r} = fscanf(fileID, '%f ', s);
				end
			end
			
			ex = fscanf(fileID, '\ndelta %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				delta = cell(1, nr);
				for r = 1:nr
					n = fscanf(fileID, '%d', 1);
					delta{r} = fscanf(fileID, '%f ', n);
				end
			end
			
			ex = fscanf(fileID, '\npsijh %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				psijh = cell(1, nr);
				for r = 1:nr
					s = fscanf(fileID, '%d %d', 2)';
					psijh{r} = fscanf(fileID, '%f ', s);
				end
			end
			
			ex = fscanf(fileID, '\nalpha %d');
			if ex
				nr = fscanf(fileID, '%d', 1);
				alpha = cell(1, nr);
				for r = 1:nr
					n = fscanf(fileID, '%d', 1);
					if n == -1
						alpha{r} = NaN;
					else
						alpha{r} = fscanf(fileID, '%d ', n);
					end
				end
			end
			
			ex = fscanf(fileID, '\nrho %d');
			if ex
				rho = fscanf(fileID, '%d', 1);
				if rho == -1
					rho = NaN;
				end
			end
			
			p = HmscPar();
			p.sigma = diag(sigma);
			p.gamma = gamma;
			p.beta = beta;
			p.V = V;
			p.nf = nf;
			p.lambda = lambda;
			p.eta = eta;
			p.delta = delta;
			p.psijh = psijh;
			p.alpha = alpha;
			p.rho = rho;
			fclose(fileID);
		end
	end
end


