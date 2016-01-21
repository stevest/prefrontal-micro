function TESTexportBackgroundStimParams(bgStimInputDend,outPath)
% Export parameter matrices in .hoc file:

% Write NMDA results to a .hoc file:
fid = fopen([outPath,'/RAPIDimportBackgroundStimParams.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref BG_Stim_Rapid[%d][%d]\n',...
    size(bgStimInputDend,1),size(bgStimInputDend,2));
fprintf(fid,'BG_rapid = %d\n',size(bgStimInputDend,2));

fprintf(fid,'\n\n// Import parameters: \n\n');
% Only for Pyramidals:
for c=1:size(bgStimInputDend,1)
    for i=1:size(bgStimInputDend,2)
        fprintf(fid,'BG_Stim_Rapid[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputDend{c,i}));
        for j=1:length(bgStimInputDend{c,i})
            fprintf(fid,'BG_Stim_Rapid[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputDend{c,i}(j)));
        end
    end
end


fclose(fid);
end