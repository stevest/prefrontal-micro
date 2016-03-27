function exportStimulationParameters(run, pathto)

% Export stimulation parameters in .hoc file:
filename = sprintf('importStimulationParameters_SN%d.hoc',run.sn);

fprintf('Exporting %s...\n',filename);
fid = fopen([pathto,filesep,filename],'W');

fprintf(fid,'// This HOC file was generated with MATLAB\n\n');
fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref PcellStimListSTR[%d]\n',run.nClusters_str);
fprintf(fid,'objref PcellStimListRND[%d]\n',run.nClusters_rnd);

fprintf(fid,'\n\n// Import parameters:\n\n');
% Structured network stimulation:
for i=1:run.nClusters_str
    fprintf(fid,sprintf('PcellStimListSTR[%d]=new Vector(%d)\n',i-1,run.nCellsInCluster_str(i)));
    for j=1:run.nCellsInCluster_str(i)
        fprintf(fid,'PcellStimListSTR[%d].x[%d]=%d\n',i-1, j-1, run.stimulatedCells_str{i}(j)-1); % NEURON idx
    end
end
% Random network stimulation:
for i=1:run.nClusters_rnd
    fprintf(fid,sprintf('PcellStimListRND[%d]=new Vector(%d)\n',i-1,run.nCellsInCluster_rnd(i)));
    for j=1:run.nCellsInCluster_rnd(i)
        fprintf(fid,'PcellStimListRND[%d].x[%d]=%d\n',i-1, j-1, run.stimulatedCells_rnd{i}(j)-1); % NEURON idx
    end
end
fprintf(fid,'//EOF\n');
fclose(fid);
end