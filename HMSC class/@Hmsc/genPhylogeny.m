function genPhylogeny(m)
if m.phylogeny == false
	error('HMSC: Model was defined without phylogeny');
end

C=eye(m.ns);
for i = 1:2:m.ns-1
	C(i,i+1)=0.99;
	C(i+1,i)=0.99;
end
m.setPhylogeny(C);

end