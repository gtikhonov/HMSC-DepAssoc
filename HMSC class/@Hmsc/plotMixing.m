function plotMixing(m, values, valuesT, dispTrue, label, showToScreen, saveToFile)

index = m.postSamInd;
if ~showToScreen
	set(0,'DefaultFigureVisible','off');
end
figure;
plot(index, values);
lim=axis;
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
xlabel('mcmc round');
tmp=ylabel(label);
set(tmp,'Interpreter','none');
tmp=title(strcat('mcmc trace plot for',{' '},label));
set(tmp,'Interpreter','none');

if(dispTrue)
	plot(index,(repmat(valuesT(:),1,length(index))),'--');
	mi=min(valuesT(:));
	ma=max(valuesT(:));
	de=(ma-mi);
	axis([lim(1),lim(2),min(lim(3),min(0,mi)-0.2*de),max(lim(4),max(0,ma)+0.2*de)]);
end

if saveToFile
	folder = strcat(m.folder,filesep,'results',filesep);
	if (~isequal(exist(folder, 'dir'),7))
		mkdir(folder);
	end
	if (~isequal(exist(folder, 'dir'),7))
		mkdir(folder);
	end
	label = strcat('mixing_', label);
	print(strcat(folder,label),'-dtiff');
end
	
if ~showToScreen
	set(0,'DefaultFigureVisible','on');
end


end
