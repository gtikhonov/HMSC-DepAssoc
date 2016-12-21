function stfMCMCRun(m, parVec, run)

posteriorFolder=strcat(m.folder,filesep,'posteriors',filesep);
[s, mess, messid] = mkdir(posteriorFolder);
fileName = strcat(posteriorFolder, 'postPar_run_', int2str(run), '.mat' );
save( fileName, 'parVec' );

end