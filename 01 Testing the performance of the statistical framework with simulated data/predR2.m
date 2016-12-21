clearvars;
globalTimer = tic();
homeFolder = '.';

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
predN = 100;
testN = 800;

NY = 100*2.^(0:5);
globalCountInt = (1:30);

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
   nyTimer = tic();
   for j = 1:length(NY)
      ny = NY(j);
      fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
      if exist(fullfile(baseFolder, 'readyPred', sprintf('readyPred-%.5d-%.3d-%.4d-%.4d.mat', ny, gC, testN, predN)), 'file')~=2
         rng(gC*10000+ny);
         load(fullfile(folder, sprintf('fm %.5d.mat', ny)), 'm')
         pred = m.predict(predN, XTest, piCellTest, [], {XrCell}, true);
         YPred = mean(cat(3,pred{:}), 3);
         save(fullfile(folder, sprintf('pred %.5d-%.4d-%.4d.mat', ny, testN, predN)), 'YPred');
         [~,~,~] = mkdir(fullfile(baseFolder, 'readyPred'));
         save(fullfile(baseFolder, 'readyPred', sprintf('readyPred-%.5d-%.3d-%.4d-%.4d.mat', ny, gC, testN, predN)), 'gC');  
      end
      toc(nyTimer);
   end
   fprintf('%s, gC = %d, true\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC);
   if exist(fullfile(baseFolder, 'readyPred', sprintf('readyPred-true-%.3d.mat', gC)), 'file')~=2
      load(fullfile(folder, sprintf('fm %.5d.mat', 100)), 'm')
      YPredTrue = nan(testN, m.ns);
      L = XTest*m.truePar.beta;
      lam = m.truePar.lambda{1};
      for i=1:testN
         lamX = lam(:,:,1) + XTest(i,2)*lam(:,:,2);
         om = lamX'*lamX;
         for j=1:m.ns
            YPredTrue(i,j) = normcdf(L(i,j), 0, sqrt(om(j,j)+1));
         end
      end
      save(fullfile(folder, sprintf('pred true.mat')), 'YPredTrue');
      [~,~,~] = mkdir(fullfile(baseFolder, 'readyPred'));
      save(fullfile(baseFolder, 'readyPred', sprintf('readyPred-true-%.3d.mat',gC)), 'gC');
   end
end
toc(globalTimer)