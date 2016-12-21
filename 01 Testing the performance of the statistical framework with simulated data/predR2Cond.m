clearvars;
addpath('C:\Google Drive\HMSC\HMSC\HMSC class');
addpath('Z:\HMSC\HMSC\HMSC class')
globalTimer = tic();
homeFolder = 'C:\DATA\HMSC\MEE revision';
% homeFolder = '.';

%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% commType = 'rich';
commType = 'poor';

% assocType = 'const';
% assocType = 'center';
% assocType = 'stress';
assocType = 'random';

% nf = 2;
nf = 5;

% prior = 'uninformed';
prior = 'informed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predCondN = 100;
predCondMCMC = 100;
testN = 800;

NY = 100*2.^(0:5);
% globalCountInt = 30:-1:1;
globalCountInt = 1:1:30;
% globalCountInt = 2+(0:15)*2;

baseFolder = fullfile(homeFolder, sprintf('%s %s %d %s', commType, assocType, nf, prior));
dataFolder = fullfile(baseFolder, 'data');
for gC = globalCountInt
   fprintf('Global counter = %d\n', gC);
   folder = fullfile(baseFolder, sprintf('gc %.3d', gC));
   load(fullfile(dataFolder, sprintf('data %.3d', gC)))
   NYMax = size(Y, 1);
   ind = mod(1:NYMax, NYMax/testN)==1;
   XTest = X(ind,:);
   YTest = Y(ind,:);
   piCellTest = piCell(ind,:);
   XrCell = [piCell, num2cell(X)];
   [~, indMax] = sort(sum(YTest), 'descend');
   indMax = indMax(1:10);
   Yc = nan(size(YTest));
   Yc(:,indMax) = YTest(:,indMax);
   nyTimer = tic();
   for j = 1:length(NY)
      ny = NY(j);
      fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
      if exist(fullfile(baseFolder, 'readyPredCond', sprintf('readyPredCond-%.5d-%.3d-%.4d-%.4d-%.4d.mat', ny, gC, testN, predCondN, predCondMCMC)), 'file')~=2
         rng(gC*10000+ny);
         load(fullfile(folder, sprintf('fm %.5d.mat', ny)), 'm')
         predCond = m.predictConditional(predCondN, Yc, predCondMCMC, XTest, piCellTest, [], {XrCell}, true);
         YPredCond = mean(cat(3,predCond{:}), 3);
         save(fullfile(folder, sprintf('condPred %.5d-%.4d-%.4d-%.4d.mat', ny, testN, predCondN, predCondMCMC)), 'YPredCond');
         [~,~,~] = mkdir(fullfile(baseFolder, 'readyPredCond'));
         save(fullfile(baseFolder, 'readyPredCond', sprintf('readyPredCond-%.5d-%.3d-%.4d-%.4d-%.4d.mat', ny, gC, testN, predCondN, predCondMCMC)), 'gC');  
      end
      toc(nyTimer);
   end
   fprintf('%s, gC = %d, true\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC);
   if exist(fullfile(baseFolder, 'readyPredCond', sprintf('readyPredCond-true-%.3d-%.4d-%.4d-%.4d.mat', gC, testN, predCondN, predCondMCMC)), 'file')~=2
      load(fullfile(folder, sprintf('fm %.5d.mat', 100)), 'm')
      rng(gC*10000);
      for i = 1:m.postSamN
         m.postSamVec(i).beta = m.truePar.beta;
         m.postSamVec(i).lambda = m.truePar.lambda;
         m.postSamVec(i).nf = m.truePar.nf;
         eta = m.postSamVec(i).eta{1};
         m.postSamVec(i).eta = {eta(:,1:m.truePar.nf)};
      end
      predCond = m.predictConditional(predCondN, Yc, predCondMCMC, XTest, piCellTest, [], {XrCell}, true);
      YPredCondTrue = mean(cat(3,predCond{:}), 3);
      save(fullfile(folder, sprintf('condPred true-%.4d-%.4d-%.4d.mat', testN, predCondN, predCondMCMC)), 'YPredCondTrue');
      [~,~,~] = mkdir(fullfile(baseFolder, 'readyPredCond'));
      save(fullfile(baseFolder, 'readyPredCond', sprintf('readyPredCond-true-%.3d-%.4d-%.4d-%.4d.mat', gC, testN, predCondN, predCondMCMC)), 'gC');
   end
end
toc(globalTimer)