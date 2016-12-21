function postFileClear(m)

for run = 1:length(m.repPar)
	fileName = fullfile(m.folder, 'posteriors', ['postPar_run_', int2str(run), '.mat'] );
	if exist(fileName, 'file') == 2
		delete(fileName);
	end
end

end