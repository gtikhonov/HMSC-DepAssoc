function z =  update_z(z, X,Xr,Xs,beta,eta,lambda,etas,lambdas,Y,nf,pi,ncr,dist,sigma,speciesX,includeXs,factorcov)
if speciesX
	[ny nc]=size(X{1});
else
	[ny nc]=size(X);
end
ns = size(Y,2);

if isempty(pi)
	nr = 0;
else
	nr = size(pi,2);
end

if speciesX
	Ez = zeros(ny,ns);
	for j=1:ns
		Ez(:,j)=X{j}*beta(:,j);
	end
else
	Ez = X*beta;
end

if includeXs == 1
	Ez = Ez + Xs*etas*lambdas;
elseif includeXs == 2
	for i=j:ns
		Ez(:,j) = Ez(:,j) + Xs{j}*etas*lambdas(:,j);
	end
end

for r=1:nr
	eta1=eta{r};
	lambda1=lambda{r};
	if factorcov(r) > 0
		Xr1 = Xr{r};
		for k=1:ncr(r)
			if factorcov(r) == 1
				XrEta = repmat(Xr1(:,k), 1, nf(r)).*eta1;
				Ez = Ez+XrEta(pi(:,r),:)*lambda1(:,:,k);
			else
				for j=1:ns
					XrEta = repmat(Xr1{j}(:,k), 1, nf(r)).*eta1;
					Ez(:,j) = Ez(:,j) + XrEta(pi(:,r),:)*lambda1(:,j,k);
				end
			end
		end
	else
		Ez = Ez+eta1(pi(:,r),:)*lambda1;
	end
end
[ny ns]=size(Ez);
if sum(dist(:,1)==3)==ns
	low=log(Y)-Ez;
	high=log(Y+1)-Ez;
elseif sum(dist(:,1)==2)==ns
	low = -Inf*ones(ny,ns);
	high = Inf*ones(ny,ns);
	cut=-Ez;
	low(Y==1)=cut(Y==1);
	high(Y==0)=cut(Y==0);
else
	low = -Inf*ones(ny,ns);
	high = Inf*ones(ny,ns);
	for j=1:ns
		if dist(j,1)==2
			cut=-Ez(:,j);
			low(Y(:,j)==1,j)=cut(Y(:,j)==1);
			high(Y(:,j)==0,j)=cut(Y(:,j)==0);
		end
		if (dist(j,1)==3)
			low(:,j)=log(Y(:,j))-Ez(:,j);
			high(:,j)=log(Y(:,j)+1)-Ez(:,j);
		end
		if (dist(j,1)==4)
			low(:,j)=Y(:,j)-Ez(:,j);
			low(Y(:,j)==0,j)=-Inf;
			high(:,j)=Y(:,j)+1-Ez(:,j);
		end
	end
end
missing=isnan(Y);
low(missing)=-Inf;
high(missing)=Inf;

%DIST=4
sel=(dist(:,1)==4);
if any(sel)
	dsigma2=diag(sigma);
	dsi2=dsigma2(sel);
	si2 = repmat(dsi2,1,ny)';
	Ez1 = Ez(:,sel);
	Y1=Y(:,sel);
	
	persistent hmc_poisson;
	if isempty(hmc_poisson)
		hmc_poisson = poisson_normal_tobit();
		hmc_poisson = hmc_poisson.input_Ey_sigma2(Y1, Ez1, si2(:,sel));
		hmc_poisson = hmc_poisson.run(1,true);
	end
	
	hmc_poisson = hmc_poisson.input_Ey_sigma2(Y1, Ez1, si2(:,sel));
	hmc_poisson = hmc_poisson.run(1,true);
	shape = size(z(:,sel));
	z(:,sel)=reshape( hmc_poisson.eta, shape);
end

% mean(hmc_poisson.eta)
%DIST=3
sel=(dist(:,1)==3);%|(dist(:,1)==3)|(dist(:,1)==4);
if any(sel)
	dsigma2=diag(sigma);
	dsi2=dsigma2(sel);
	si2 = repmat(dsi2,1,ny)';
    
    si2( si2>1E2)=1E2;
	si2(si2<1E-3)=1E-3;

	
	Ez1 = Ez(:,sel);
	Y1=Y(:,sel);

    r = 1000;
%     r = exp(Ez1)*100;
%     r(r>1) = 1;
	logr = log(r);
  	
	missing=isnan(Y1);
  	Y1(missing) = poissrnd(exp(Ez1(missing)));

  %posterior w
	Psi = z(:,sel) - logr;
	
	n1 = size(Y1,1);
	n2 = size(Y1,2);
	
	b =  r;
	c =  abs(Psi);
	m= b ./2 ./c .* tanh(c/2);
	
	%when c is too big, this converges to 2;
	v_part2 =   (sinh(c)-c) ./ cosh(c ./2).^2;
	v_part2(isnan(v_part2))=2;
	
	v= b ./4 ./c.^3 .* v_part2;
	
	w = randn(n1,n2) .* sqrt(v) + m;
    
    w = abs(w);
    
	% 	w(w<0) =0;
	
	%posterior Psi
   
	Psi_var = 1.0 ./ (w + 1./si2);
	
	
	Psi_mean = Psi_var .* ( Y1 - r./2.0 + (Ez1 - logr) ./ si2);
	Psi = normrnd( Psi_mean, sqrt(Psi_var));
	z(:,sel) =  Psi + logr;
end
%DIST=2
sel=(dist(:,1)==2);
if any(sel)
	mu = zeros(ny,sum(sel));
	dsigma=sqrt(diag(sigma));
	dsi=dsigma(sel);
	si = repmat(dsi,1,ny)';
	eps = Hmsc.tnormrnd(mu,si,low(:,sel),high(:,sel));
	[row,col]=find(eps==Inf);
	for  j=1:length(row)
		eps(row(j),col(j))=low(row(j),col(j))+0.5;
	end
	[row,col]=find(eps==-Inf);
	for  j=1:length(row)
		eps(row(j),col(j))=high(row(j),col(j))-0.5;
	end
	z(:,sel)=Ez(:,sel)+eps;
end
% IF NUMERICAL ERRORS APPEAR, USE THIS VERSION FOR DIST NOT 4
% z(:,sel)=Ez(:,sel)+max(min(eps,10),-10);

sel=(dist(:,1)==1);
if any(sel)
	dsigma=sqrt(diag(sigma));
	Y1=Y(:,sel);
	missing=isnan(Y1);
	Ez1=Ez(:,sel);
	dsi=dsigma(sel);
	mu = zeros(ny,sum(sel));
	si = repmat(dsi,1,ny)';
	low = -Inf(ny,sum(sel));
	high = Inf(ny,sum(sel));
	eps = Hmsc.tnormrnd(mu(missing),si(missing),low(missing),high(missing));
	z1 = Y1;
	z1(missing) = Ez1(missing)+eps;
	z(:,sel) = z1;
	
end

