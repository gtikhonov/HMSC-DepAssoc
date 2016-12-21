clearvars;
repN = 1000;

calcFlag = false;

folder = '.';

file=strcat('.','\data\data Y.xlsx');
da=importdata(file);
Y = da.data.theData;
spNames = da.colheaders.theData;
n = size(Y, 1);

file=strcat('.','\data\data X.xlsx');
da=importdata(file);
X = da.data.theData;
X = [ones(n, 1), X];
X(:,3) = log(X(:,3)+1);
X(:,4) = X(:,3).*X(:,2);
covNames =  ['intercept', da.colheaders.theData, 'SmDist'];

file=strcat('.','\data\data Pi.xlsx');
da=importdata(file);
Pi = da.data;
piCell = arrayfun(@(c)sprintf('%.3d', c), Pi, 'UniformOutput', false);

m = Hmsc(folder, false, false, false, 0, false, [false,false], [true,false]);

XrCell = [piCell(:,1),num2cell(X(:,1:3))];
m.setData(Y,'probit',X,piCell,[],[],[],{XrCell,[]},[]);
m.setCovNames(covNames);
m.setSpeciesNames(spNames);
m.setCovScaling([2,1,1,1]);
m.setFactorCovScaling({[2,1,1]});
m.setPriorsDefault();
m.setMCMCOptions(1000, 10);
m.setMCMCAdapt([repN*0.3,0], false, [0,0], false);
m.setMCMCSaveOptions(true, false);
m.sampleMCMC(repN, false, [], 3);
m.setPostThinning((repN*0.4):m.repN,5);
m.postRamClear();

