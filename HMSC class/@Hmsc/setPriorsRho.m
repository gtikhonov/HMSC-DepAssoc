function setPriorsRho(m, rhop, rhow)

if isempty(rhop) && ~isempty(rhow)
	error('HMSC: Descrete grid must be explicitely specified when explicitely setting the weights');
end

if ~isempty(rhop)
	if any(rhop<-1) || any(rhop>1)
		error('HMSC: Descrete grid must be within [-1; 1] interval');
	end
else
	rhop =  -1:0.01:1;
end
if ~isempty(rhow)
	if length(rhow) ~= length(rhop)
		error('HMSC: Length of weights vector must be equal to the length of grid vector');
	end
	if any(rhow) < 0
		error('HMSC: Weights must be positive');
	end
	if abs(sum(rhow)-1) > 0.00001
		error('HMSC: Sum of weights must be equal to 1');
	end
else
	rhow = ones(1,length(rhop));		
	ind0 = rhop==0;
	if any(ind0)
		rhow = 0.5*rhow/(sum(rhow)-1);
		rhow(ind0) = 0.5;
	else
		rhow = rhow/sum(rhow);
	end
end
m.rhopw = [rhop; rhow]';

end