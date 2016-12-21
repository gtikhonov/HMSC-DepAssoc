clearvars;
homeFolder = '.';
globalTimer = tic();
ns = 50;

%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% commType = 'rich';
commType = 'poor';

% assocType = 'const';
% assocType = 'center';
assocType = 'stress';
% assocType = 'random';

nf = 2;
% nf = 5;

% prior = 'uninformed';
prior = 'informed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NY = 100*2.^(0:5);
NYMax = 100*2.^(7);
nyN = length(NY);

globalCountInt = (1:30);

baseFolder = fullfile(homeFolder, sprintf('%s %s %d %s', commType, assocType, nf, prior));
minX = -1;
maxX = 1;
cenX = (minX+maxX)/2;
threshold = 0.95;
trilVec = tril(ones(ns), -1);
trilVec = (trilVec(:)>0);


if exist(fullfile(baseFolder, 'infSumSDM.mat'), 'file')==2
   load(fullfile(baseFolder, 'infSum.mat'), 'matShift');
   load(fullfile(baseFolder, 'infSumSDM.mat'), 'matShiftSDM')
else
   matShift = nan(nyN, 2, length(globalCountInt));
   matShiftSDM = matShift;
end

for gC = globalCountInt
   fprintf('Global counter = %d\n', gC);
   folder = fullfile(baseFolder, sprintf('gc %.3d', gC));
   if exist(fullfile(baseFolder, 'infSumSDM.mat'), 'file')==2
      matRepShift = permute(matShift(:,gC,:), [1,3,2]);
      matRepShiftSDM = permute(matShiftSDM(:,gC,:), [1,3,2]);
      save(fullfile(folder, 'infSumRep.mat'), 'matRepShift');
      save(fullfile(folder, 'infSumRepSDM.mat'), 'matRepShiftSDM');
      [~,~,~] = mkdir(fullfile(baseFolder, 'readyInf'));
      save(fullfile(baseFolder, 'readyInf', sprintf('readyInf-%.3d.mat', gC)), 'gC');
      continue;
   end
   
   if exist(fullfile(baseFolder, 'readyInf', sprintf('readyInf-%.3d.mat', gC)), 'file')~=2
      matRepShift = nan(nyN, 2);
      matRepShiftSDM = matRepShift;
      for j = 1:length(NY)
         nyTimer = tic();
         ny = NY(j);
         fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
         rng(gC*10000+ny);
         load(fullfile(folder, sprintf('fm %.5d.mat', ny)), 'm')
         low = m.getPostOmega(1,[1,minX]);
         lowT = low{2};
         low = low{1};
         cen = m.getPostOmega(1,[1,cenX]);
         cenT = cen{2};
         cen = cen{1};
         high = m.getPostOmega(1,[1,maxX]);
         highT = high{2};
         high = high{1};
         low = cell2mat(cellfun(@(c) reshape(corrcov(reshape(c, m.ns, m.ns)), 1, m.ns^2), num2cell(low, 2), 'UniformOutput' , false));
         cen = cell2mat(cellfun(@(c) reshape(corrcov(reshape(c, m.ns, m.ns)), 1, m.ns^2), num2cell(cen, 2), 'UniformOutput' , false));
         high = cell2mat(cellfun(@(c) reshape(corrcov(reshape(c, m.ns, m.ns)), 1, m.ns^2), num2cell(high, 2), 'UniformOutput' , false));
         lowT = corrcov(reshape(lowT, m.ns, m.ns));
         cenT = corrcov(reshape(cenT, m.ns, m.ns));
         highT = corrcov(reshape(highT, m.ns, m.ns));
         cmpT = sign(highT-lowT);
         
         warning('off','all')
         sSDM = nan(m.ns);
         for sp1 = 1:m.ns
            for sp2 = 1:m.ns
               if sp1~= sp2
                  sdmX = [m.X(:,2), m.X.*repmat(m.Y(:,sp2),1,2)];
                  [b,~,stats] = glmfit(sdmX, m.Y(:,sp1),  'binomial', 'Link', 'probit');
                  sSDM(sp1,sp2) = (stats.p(4) < 0.05)*sign(b(4));
               end
            end
         end
         warning('on','all')
         estTrue = sum(sum((sSDM==cmpT)));
         comp = cmpT&(sSDM==-cmpT)|(~cmpT)&(sSDM~=0);
         comp = comp - diag(diag(comp));
         estWrong = sum(sum(comp));
         matRepShiftSDM(j,:) = [estTrue, estWrong];
         
         cmp = reshape( sum(low < high)/m.postSamN, m.ns, m.ns);
         cmp = -(cmp<(1-threshold)) + (cmp>threshold);
         cmp = cmp(trilVec);
         cmpT = cmpT(trilVec);
         estTrue = sum((cmp==cmpT));
         estWrong = sum(cmpT&(cmp==-cmpT)|(~cmpT)&cmp);
         matRepShift(j, :) = [estTrue, estWrong];
         
         toc(nyTimer);
      end
      save(fullfile(folder, 'infSumRep.mat'), 'matRepShift');
      save(fullfile(folder, 'infSumRepSDM.mat'), 'matRepShiftSDM');
      [~,~,~] = mkdir(fullfile(baseFolder, 'readyInf'));
      save(fullfile(baseFolder, 'readyInf', sprintf('readyInf-%.3d.mat', gC)), 'gC');
   else
      load(fullfile(folder, 'infSumRep.mat'), 'matRepShift');
      load(fullfile(folder, 'infSumRepSDM.mat'), 'matRepShiftSDM');
   end
   matShift(:,:,gC) = matRepShift;
   matShiftSDM(:,:,gC) = matRepShiftSDM;
   toc(globalTimer)
end


mat = matShift;
figure;
tabShift = (nanmean(mat, 3));
tabShiftSd = (nanstd(mat, 0, 3));
tabShift = [tabShift, (ns^2-ns)/2 - sum(tabShift, 2)] / ((ns^2-ns)/2)*100;
tabShiftSd = [tabShiftSd, zeros(nyN,1)] / ((ns^2-ns)/2)*100;
b = bar(1:nyN,tabShift);
hold on;
b(1).FaceColor = 'red';
b(2).FaceColor = 'blue';
b(3).FaceColor = [0.8, 0.8, 0.8];
errorbar((1:nyN)-0.22,tabShift(:,1),tabShiftSd(:,1), 'k.', 'color', 'black', 'LineWidth', 2 )
errorbar((1:nyN),tabShift(:,2),tabShiftSd(:,2),  'k.', 'color', 'black', 'LineWidth', 2 )
ylim([0,105])
xlim([0.5,6.5])
set(gca,'FontSize',18)
set(gca,'Xtick', [1:6], 'XTickLabel', NY);

mat = matShiftSDM;
figure;
tabShift = (nanmean(mat, 3));
tabShiftSd = (nanstd(mat, 0, 3));
tabShift = [tabShift, (ns^2-ns) - sum(tabShift, 2)] / ((ns^2-ns))*100;
tabShiftSd = [tabShiftSd, zeros(nyN,1)] / ((ns^2-ns))*100;
b = bar(1:nyN,tabShift);
hold on;
b(1).FaceColor = 'red';
b(2).FaceColor = 'blue';
b(3).FaceColor = [0.8, 0.8, 0.8];
errorbar((1:nyN)-0.22,tabShift(:,1),tabShiftSd(:,1), 'k.', 'color', 'black', 'LineWidth', 2 )
errorbar((1:nyN),tabShift(:,2),tabShiftSd(:,2),  'k.', 'color', 'black', 'LineWidth', 2 )
ylim([0,105])
xlim([0.5,6.5])
set(gca,'FontSize',18)
set(gca,'Xtick', [1:6], 'XTickLabel', NY);


