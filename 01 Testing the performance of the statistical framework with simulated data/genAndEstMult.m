clearvars;
globalTimer = tic();

%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
commType = 'rich';
% commType = 'poor';

% assocType = 'const'
% assocType = 'center';
assocType = 'stress';
% assocType = 'random';

% nf = 2;
nf = 5;

% prior = 'uninformed';
% prior = 'informed';
prior = 'correct';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcThin = 10;
mcmcRep = 10;
mcmcSam = 100;
mcmcBurn = 5;
mcmcPostBurn = 5;
mcmcPostThin = 5;

% NY = 100*2.^(0:1);
NY = 100*2.^(0:5);
NYMax = 100*2.^(7);

% globalCountInt = 1*60+(1:60);
globalCountInt = (1:15);

for gC = globalCountInt
   fprintf('Global counter = %d\n', gC);
   baseFolder = sprintf('%s %s %d %s', commType, assocType, nf, prior);
   [~,~,~] = mkdir(baseFolder);
   dataFolder = fullfile(baseFolder, 'data');
   [~,~,~] = mkdir(dataFolder);
   folder = fullfile(baseFolder, sprintf('gc %.3d', gC));
   genBetaLambda(gC, commType, assocType, nf, dataFolder)
   gen(NYMax, nf, gC, dataFolder);
   
   for j = 1:length(NY)
      nyTimer = tic();
      ny = NY(j);
      fprintf('%s, gC = %d, ny = %d\n', sprintf('%s %s %d %s', commType, assocType, nf, prior), gC, ny);
      if exist(fullfile(baseFolder, 'ready', sprintf('ready-%.5d-%.3d.mat', ny, gC)), 'file')~=2
         est(ny, NYMax, gC, nf, prior, mcmcSam, mcmcThin, mcmcBurn, mcmcRep, mcmcPostBurn, mcmcPostThin, folder, dataFolder);
         saveReady(ny, gC, baseFolder);
      end
      toc(nyTimer);
   end
end
toc(globalTimer)
sendEmail