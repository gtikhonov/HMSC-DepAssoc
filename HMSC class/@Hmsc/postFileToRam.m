function postFileToRam(m, loadVec)

for run=loadVec
	fileName = fullfile(m.folder, 'posteriors', ['postPar_run_', int2str(run), '.mat'] );
	load( fileName, 'parVec' );
	m.repPar{run} = parVec;
end

end