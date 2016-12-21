function plotBox(m, values, valuesT, dispTrue, label, showToScreen, saveToFile)
if ~showToScreen
	set(0,'DefaultFigureVisible','off');
end

figure;
boxplot(values);
lim=axis;
hold on;
tmp=ylabel(label);
set(tmp,'Interpreter','none');
if(dispTrue)
	scatter(1:length(valuesT),valuesT,'filled','red');
	mi=min(valuesT(:));
	ma=max(valuesT(:));
	de=(ma-mi);
	axis([lim(1),lim(2),min(lim(3),min(0,mi)-0.2*de-0.2),max(lim(4),max(0,ma)+0.2*de)+0.2]);
end

if saveToFile
	folder = strcat(m.folder,filesep,'results',filesep);
	if (~isequal(exist(folder, 'dir'),7))
		mkdir(folder);
	end
	if (~isequal(exist(folder, 'dir'),7))
		mkdir(folder);
	end
	label = strcat('box_', label);
	print(strcat(folder,label),'-dtiff');
end
	
if ~showToScreen
	set(0,'DefaultFigureVisible','on');
end



% nboxes=nrow*ncol;
% cols=ceil(sqrt(nboxes)/ncol);
% rows=ceil(nrow/cols);
% for i=1:nrow
% 	h=subplot(rows,cols,i);
% 	p=get(h,'pos');
% 	boxplot(values(:,(i-1)*ncol+1:(i-1)*ncol+ncol));
% 	lim=axis;
% 	hold on;
% 	tmp=ylabel(lab);
% 	set(tmp,'Interpreter','none');
% 	%set(gca,'xtick',1:ncol,'xticklabel',xticks);
% 	set(gca,'xtick',1:ncol);
% 	% set(gca,'xticklabel',xticks);
% 	tmp=title(strcat(ylab,yticks{i}));
% 	set(tmp,'Interpreter','none');
% 	if(i==nrow)
% 		tmp=xlabel(xlab);
% 		set(tmp,'Interpreter','none');
% 	end
% 	plot([0.5 ncol+0.5],[0 0],'blue');
% 	if true_values
% 		scatter(1:length(valuesT(:,i)),valuesT(:,i),'filled','red');
% 		mi=min(valuesT(:,i));
% 		ma=max(valuesT(:,i));
% 		de=(ma-mi);
% 		axis([0.5,ncol+0.5,min(lim(3),min(0,mi)-0.2*de-0.2),max(lim(4),max(0,ma)+0.2*de)+0.2]);
% 	end
% end
