function est(ny, NYMax, globalCount, nf, prior, mcmcSam, mcmcThin, mcmcBurn, mcmcRep, mcmcPostBurn, mcmcPostThin, folder, dataFolder)
rng(globalCount);
load( fullfile(dataFolder, sprintf('data %.3d.mat', globalCount) ))

ind = mod(1:NYMax, NYMax/ny)==0;
X = X(ind,:);
Y = Y(ind,:);
piCell = piCell(ind,:);

m = Hmsc(folder, false, false, false, 0, false, [false], [1]);
m.setData(Y, 'probit', X, piCell, [], [], [], [], [])
XrCell = [piCell, num2cell(X)];
m.setRanFactorCov({XrCell});
m.setFactorCovScaling({[0,0]});
m.setPriorsDefault()


m.setMCMCOptions(mcmcSam, mcmcThin);
m.setMCMCSaveOptions(true, false);
switch prior
   case 'uninformed'
      m.setMCMCAdapt([mcmcBurn, 0], false);
   case 'informed'
      m.setMCMCAdapt([mcmcBurn, 0], nf);
      nur = 3;
      a1r = 1; b1r = 1;
      a2r = 1; b2r = 1;
      m.setPriorsRandomFactors(nur, a1r, b1r, a2r, b2r)
   otherwise
      error('Wrong prior type')
end
m.sampleMCMC(mcmcRep, false, [], 10);
m.setPostThinning(mcmcPostBurn:m.repN, mcmcPostThin);
m.postRamClear()

m.setParameters(tP);
[~,~,~] = mkdir(folder);
save(fullfile(folder, sprintf('fm %.5d.mat', ny)), 'm')

