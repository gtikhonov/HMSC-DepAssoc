clearvars;
globalTimer = tic();
homeFolder = '.';
ns = 50;

%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% commType = 'rich';
commType = 'poor';

% assocType = 'const'
% assocType = 'center';
assocType = 'stress';
% assocType = 'random';


nf = 2;
% nf = 5;

% prior = 'uninformed';
prior = 'informed';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testN = 800;

NY = 100*2.^(0:5);
globalCountInt = (15:30);

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
   
   nyTimer = tic();
   warning('off','all')
   for j = 1:length(NY)
      ny = NY(j);
      ind = mod(1:NYMax, NYMax/ny)==0;
      XTrain = X(ind,:);
      YTrain = Y(ind,:);
      fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
      YPredSDM = nan(testN, ns);
      for j = 1:ns
         b = glmfit(XTrain(:,2), YTrain(:,j), 'binomial', 'Link', 'probit');
         YPredSDM(:,j) = glmval(b, XTest(:,2), 'probit');
      end
      save(fullfile(folder, sprintf('predSDM %.5d-%.4d.mat', ny, testN)), 'YPredSDM');
      toc(nyTimer);
   end
end
warning('on','all')
toc(globalTimer)
