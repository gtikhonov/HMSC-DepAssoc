function plotVariancePartitioning(m,fixed, fixedsplit,random, traitR2,groupnames,showToScreen, saveToFile);
if ~showToScreen
	set(0,'DefaultFigureVisible','off');
end
ngroups=size(fixedsplit,2);
fixedsplit2=repmat(fixed,1,ngroups).*fixedsplit;
avars=[fixedsplit2 random];
%[tmp,index]=sort(fixex);
%index=flip(index);
index=1:m.ns;
figure;
if m.ns>1
	bar(avars(index,:),'stacked');
else
	bar([avars; zeros(1,size(avars,2))],'stacked') ;
end
xlim([0 m.ns+1]);
xlabel('species');
ylabel('proportion of variance');
ylim([0 1]);
title(['Variance partitioning. Traits explain ',num2str(100*traitR2,2),'% of fixed effects']);
set(gca,'fontsize',16)
	if m.ns>1
	me=mean(avars);
else
	me=avars;
end
jj=0;
for i=1:ngroups
	jj=jj+1;
	name=['fixed: ',groupnames{i}];
	legendInfo{jj}=[name,' (mean = ',num2str(100*me(jj),2),'%)'];
end
for level=1:m.nr
	jj=jj+1;
	name=['random: ',m.levelNames{level}];
	legendInfo{jj}=[name,' (mean = ',num2str(100*me(jj),2),'%)'];
end
if saveToFile
	folder = strcat(m.folder,filesep,'results',filesep);
	if (~isequal(exist(folder, 'dir'),7))
		mkdir(folder);
	end
	print(strcat(folder,'variance-partitioning.tiff'),'-dtiff');
	legend(legendInfo,'Location','southoutside','interpreter','none');
	print(strcat(folder,'variance-partitioning_legend.tiff'),'-dtiff');
end

if ~showToScreen
	set(0,'DefaultFigureVisible','on');
end