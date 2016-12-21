function setPriorsDefault(m)

m.nu = 4; % gamma parameter for outlierspecies

m.nur = 3;                                   % gamma hyperparameters for psijh
m.a1r = 50.; m.b1r = 1;                         % gamma hyperparameters for delta_ 1
m.a2r = 50.; m.b2r = 1;                         % gamma hyperparameters delta_h, h >= 2

m.nus = 3;                                   % gamma hyperparameters for psijh
m.a1s = 50.; m.b1s = 1;                         % gamma hyperparameters for delta_ 1
m.a2s = 50.; m.b2s = 1;                         % gamma hyperparameters delta_h, h >= 2

m.asigma = 1.0*ones(m.ns,1);                   % gamma hyperparameters sigma
m.bsigma = 0.3*ones(m.ns,1);

m.f0 = m.nc + 1;
m.V0 = eye(m.nc);
m.Ugamma = eye(m.nc*m.nt);
m.mgamma = zeros(m.nc*m.nt,1);

m.alphapw = cell(1, m.nr);
for r = 1:m.nr
	alphapw=[];
	if m.spatial(r)
		%d = sqrt( sum((max(m.xy{r},1)-min(m.xy{r},1)).^2) );
		d = sqrt(sum((max(m.xy{r})-min(m.xy{r})).^2));
		alphap = (0:0.01:1)*d;
		alphaw = ones(1,length(alphap));
		alphaw = 0.5*alphaw/(sum(alphaw)-1);
		alphaw(1) = 0.5;
		alphapw = [alphap; alphaw]';
	end
	m.alphapw{r} = alphapw;
end

rhop =  -1:0.01:1;
rhow = ones(1,length(rhop));
rhow = 0.5*rhow/(sum(rhow)-1);
ind0 = rhop==0;
rhow(ind0) = 0.5;
m.rhopw = [rhop; rhow]';

end
