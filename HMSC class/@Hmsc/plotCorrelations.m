function plotCorrelations(m,correlations, index, plottitle, type, showToScreen, saveToFile)

if ~showToScreen
	set(0,'DefaultFigureVisible','off');
end
if type(1)
	figure;
	imagesc(correlations(index,index));
% 	colorbar
	greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
	redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
	colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';
	colorMap = [redColorMap; zeros(1, 256); greenColorMap]';
	colorMap = [1-redColorMap; 1-redColorMap-greenColorMap; 1-greenColorMap]';
	colormap(colorMap);
	ma=max(max(abs(correlations(:))),1);
	caxis([-ma ma]);
	set(gca,'fontsize',16)
% 	title(plottitle,'interpreter','none');
	if saveToFile
		folder = strcat(m.folder,filesep,'results',filesep);
		if (~isequal(exist(folder, 'dir'),7))
			mkdir(folder);
		end
		print(strcat(folder,'matrix_',plottitle,'.tiff'),'-dtiff');
	end
end
if type(2)
	n=length(correlations);
	myLabel = cell(n);
	for i = 1:n
		myLabel{i} = num2str(index(i));
	end
	pos=ones(n,n).*(correlations>0);
	neg=ones(n,n).*(correlations<0);
	myColorMapPos = zeros(n,3);
	myColorMapNeg = zeros(n,3);
	myColorMapPos(:,1)=1;
	myColorMapNeg(:,3)=1;
	figure;
	circularGraph(neg(index,index),'Colormap',myColorMapNeg,'Label',myLabel);
	circularGraph(pos(index,index),'Colormap',myColorMapPos,'Label',myLabel);
	% title(plottitle);
	if saveToFile
		folder = strcat(m.folder,filesep,'results',filesep);
		if (~isequal(exist(folder, 'dir'),7))
			mkdir(folder);
		end
		print(strcat(folder,'circle_',plottitle,'.tiff'),'-dtiff');
	end
end


if ~showToScreen
	set(0,'DefaultFigureVisible','on');
end
