function exportNetworkParameters(run, conf, pathto)

if strcmp(conf,'rnd')
    s=0;
elseif strcmp(conf,'str')
    s=1;
else
    error('No network configuration given! Exiting.');
    return;
end
% Export Network Connectivity:
% Export STRUCTURED parameter matrices in .hoc file:

filename = sprintf('importNetworkParameters%s_SN%d.hoc',upper(conf),run.sn);

fprintf('Exporting %s...\n',filename);
fid = fopen([pathto,filesep,filename],'W');

fprintf(fid,'// This HOC file was generated with MATLAB\n');
fprintf(fid,sprintf('nPCcells=%d\n',run.nPC));
fprintf(fid,sprintf('nPVcells=%d\n',run.nPV));
fprintf(fid,sprintf('nAllCells=%d\n',run.nAll));
fprintf(fid,sprintf('steps_per_ms=%d\n',run.stepsperms));

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref C, W\n');
fprintf(fid,'C = new Matrix(nAllCells, nAllCells)\n');
fprintf(fid,'W = new Matrix(nPCcells, nPCcells)\n');

fprintf(fid,'\n// Import parameters: (long-long text following!)\n');
% network connectivity:
for i=1:run.nAll
    for j=1:run.nAll
        if s
            tmp = run.configuration_str(i,j);
        else
            tmp = run.configuration_rnd(i,j);
        end
        fprintf(fid,'C.x[%d][%d]=%d\n',i-1,j-1, tmp);
    end
end
% Network synaptic weights
for i=1:run.nPC
    for j=1:run.nPC
        if s
            tmp = run.weights_str(i,j);
        else
            tmp = run.weights_rnd(i,j);
        end
        fprintf(fid,'W.x[%d][%d]=%f\n',i-1,j-1, tmp);
    end
end
fprintf(fid,'//EOF\n');
fclose(fid);

end