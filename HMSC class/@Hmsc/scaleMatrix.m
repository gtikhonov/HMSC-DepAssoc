function [scale, XSc] = scaleMatrix(X, scaleFlag, speciesX)
%SCALEMATRIX Summary of this function goes here
%   Detailed explanation goes here
if speciesX == false
	nc = size(X,2);
else
	ns = length(X);
	nc = size(X{1},2);
end
if( length(scaleFlag)~=nc )
	error('HMSC: length of scaling flags vector must be equal to number of covariates');
end
ind = scaleFlag==2;
if any(scaleFlag < 0) || any(scaleFlag > 2) || sum(ind)>1
	error('HMSC: elements of scaling flags could take integer values from 0 to 2 and not more than one must be marked as intercept');
end

scale = NaN(2, nc);
XSc = X;
if any(scaleFlag~=0)
	if speciesX == false
		XAll = X;
	else
		XAll = cat(1, X{1:ns});
	end
	for i = 1:nc
		if scaleFlag(i) == 1
			if sum(ind) == 1 % intercept specified, centering and scaling
				me = mean(XAll(:,i));
			else % no intercept, only scaling
				me = 0;
			end
			sd = sqrt(sum((XAll(:,i)-me).^2)/(size(XAll,1)-1));
			scale(1:2,i) = [me; sd];
			if speciesX == false
				XSc(:,i) = (X(:,i)-me)/sd;
			else
				for j = 1:ns
					XSc{j}(:,i) = (X{j}(:,i)-me)/sd;
				end
			end
		end
	end
end

end

