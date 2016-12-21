function postRamToFile(m, loadVec)

[~, ~, ~] = mkdir(posteriorFolder);
for run=loadVec
	fileName = fullfile(m.folder, 'posteriors', ['postPar_run_', int2str(run), '.mat'] );
	parVec = m.repPar{run};
	save( fileName, 'parVec' );
end

end