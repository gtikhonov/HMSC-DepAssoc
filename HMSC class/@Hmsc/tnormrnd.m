function eps = tnormrnd(mu,si,low,high)

[ny, ns]=size(mu);
tmp = unifrnd(0,1,[ny ns]);
eps=mu+sqrt(2).*si.*erfinv(-((-1+tmp).*erf((low-mu)./(sqrt(2)*si)))+tmp.*erf((-mu+high)./(sqrt(2)*si)));

end