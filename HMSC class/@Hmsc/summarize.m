function res=summarize(m, name, values, valuesT, q, dim1, dim2, tran, rNames, cNames, dispTrue, showToScreen, saveToFile)

res = cell(1, length(q)+dispTrue+1);
for k=0:(length(q)+dispTrue)
	if k==0
		vals=reshape(mean(values),dim1,dim2);
	elseif k<=length(q)
		vals=reshape(quantile(values,q(k)),dim1,dim2);
	else
		vals=reshape(valuesT,dim1,dim2);
	end
	if tran
		[d2, d1] = deal(dim1, dim2);
		vals = transpose(vals);
	else
		[d1, d2] = deal(dim1, dim2);
	end
	tab = cell2table(num2cell(vals));
	if ~isempty(cNames)
		tab.Properties.VariableNames = cNames;
	end
	if ~isempty(rNames)
		tab.Properties.RowNames = rNames;
	end
	
	if saveToFile
		if k == 0
			fName = sprintf('%s mean.csv', name); 
		elseif k <= length(q)
			fName = sprintf('%s quantile %s.csv', name, num2str(q(k)));
		else
			fName = sprintf('%s true.csv', name);
		end
		fileName = fullfile(m.folder,'results',fName);
		writetable(tab, fileName, 'WriteRowNames',true);
	end	
	
	if showToScreen
		if k == 0
			fprintf('Mean of %s\n', name);
		elseif k <= length(q)
			fprintf('Quantile %s of %s\n', num2str(q(k)), name);
		else
			fprintf('True values of %s\n', name);
		end
		disp(tab);
	end
	
% 		if tran
% 		[d2, d1] = deal(dim1, dim2);
% 		vals = transpose(vals);
% 	else
% 		[d1, d2] = deal(dim1, dim2);
% 	end
% 	if isempty(rNames)
% 		sh2 = 0;
% 	else
% 		sh2 = 1;
% 	end;
% 	if isempty(cNames)
% 		sh1 = 0;
% 	else
% 		sh1 = 1;
% 	end;
% 	
% 	A = cell(d1+sh1,d2+sh2);
% 	if sh1 && sh2
% 		A{1,1}='';
% 	end
% 	if sh1
% 		for j=1:d2
% 			A{1,j+sh2}=cNames{j};
% 		end
% 	end
% 	if sh2
% 		for i=1:d1
% 			A{i+sh1,1}=rNames{i};
% 		end
% 	end
% 	for i=1:d1
% 		for j=1:d2
% 			A{i+sh1,j+sh2}=num2str(vals(i,j), '%f');
% 		end
% 	end
% 	
% 	if saveToFile
% 		fileName = strcat(m.folder,filesep,'results',filesep,name,'.xlsx');
% 		if k == 0
% 			[status, message] = xlswrite(fileName,A, 1);
% 		elseif k <= length(q)
% 			[status, message] = xlswrite( fileName,A, sprintf('q %s', num2str(q(k))) );
% 		else
% 			[status, message] = xlswrite(fileName,A,'True');
% 		end
% 		if status ~= 1
% 			error( message );
% 		end
% 	end
% 	if showToScreen
% 		if k == 0
% 			fprintf('Mean of %s\n', name);
% 		elseif k <= length(q)
% 			fprintf('Quantile %s of %s\n', num2str(q(k)), name);
% 		else
% 			fprintf('True values of %s\n', name);
% 		end
% 		disp(A);
% 	end
	res{k+1} = tab;
end

