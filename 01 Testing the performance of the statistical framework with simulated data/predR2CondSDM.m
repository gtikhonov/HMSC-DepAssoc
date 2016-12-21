clearvars;
globalTimer = tic();
homeFolder = '.';
ns = 50;

%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% commType = 'rich';
commType = 'poor';

% assocType = 'const'
% assocType = 'center';
% assocType = 'stress';
assocType = 'random';

% nf = 2;
nf = 5;

% prior = 'uninformed';
prior = 'informed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testN = 800;

NY = 100*2.^(0:5);
globalCountInt = (16:30);

for gC = globalCountInt
   fprintf('Global counter = %d\n', gC);
   baseFolder = fullfile(homeFolder, sprintf('%s %s %d %s', commType, assocType, nf, prior));
   dataFolder = fullfile(baseFolder, 'data');
   folder = fullfile(baseFolder, sprintf('gc %.3d', gC));
   load(fullfile(dataFolder, sprintf('data %.3d', gC)))
   NYMax = size(Y, 1);
   ind = mod(1:NYMax, NYMax/testN)==1;
   XTest = X(ind,:);
   YTest = Y(ind,:);
   [~, indMax] = sort(sum(YTest), 'descend');
   indMax = sort(indMax(1:10));
   Yc = nan(size(YTest));
   Yc(:,indMax) = YTest(:,indMax);
   nyTimer = tic();
   
   warning('off','all')
   for j = 1:length(NY)
      ny = NY(j);
      ind = mod(1:NYMax, NYMax/ny)==0;
      XTrain = X(ind,:);
      YTrain = Y(ind,:);
      YTestSDM = YTest(:,indMax);
      sdmXTest = [XTest(:,2), YTestSDM, repmat(XTest(:,2), 1, size(YTestSDM,2)).*YTestSDM];
      sdmY = YTrain(:,indMax);
      sdmX = [XTrain(:,2), sdmY, repmat(XTrain(:,2), 1, size(sdmY,2)).*sdmY];
      fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
      YPredCondSDM = nan(testN, ns);
      for j = 1:ns
         if ~ismember(j, indMax)
            b = glmfit(sdmX, YTrain(:,j), 'binomial', 'Link', 'probit');
            YPredCondSDM(:,j) = glmval(b, sdmXTest, 'probit');
         end
      end
      YPredCondSDM(:,indMax) = YTest(:,indMax);
      save(fullfile(folder, sprintf('condPredSDM %.5d-%.4d.mat', ny, testN)), 'YPredCondSDM');
      toc(nyTimer);
   end
end
warning('on','all')
toc(globalTimer)
