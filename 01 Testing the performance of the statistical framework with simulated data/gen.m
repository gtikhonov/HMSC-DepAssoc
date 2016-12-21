function gen(ny, nf, globalCount, folder)
rng(globalCount);

ns = 50;
nc = 2;
X = [ones(ny, 1), linspace(-1+2/ny, 1, ny)'];

m = Hmsc(folder, false, false, false, 0, false, [false], [1]);
m.setDim(ny, ns, nc)
m.setX(X);
m.genPi([m.ny]);
m.setDist('probit');
XrCell = [m.piCell, num2cell(X)];
m.setRanFactorCov({XrCell});
m.setPriorsDefault();
tP = HmscPar();
tP.nf = nf;
m.setParameters(tP);
m.genParameters();
load(fullfile(folder, sprintf('betaLambda %.3d.mat', globalCount)))
tP = m.truePar;
tP.beta = betaT;
tP.lambda = {lambdaT};
m.setParameters(tP);
m.genY(0);

Y = m.Y;
piCell = m.piCell;
save(fullfile(folder, sprintf('data %.3d.mat', globalCount)), 'Y', 'X', 'piCell', 'tP')

