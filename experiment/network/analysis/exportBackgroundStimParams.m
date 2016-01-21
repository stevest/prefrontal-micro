function exportBackgroundStimParams(bgStimInputDend,bgStimInputApic,bgStimInputPV, bgStimInputCB, bgStimInputCR,outPath)
% Export parameter matrices in .hoc file:

% Write NMDA results to a .hoc file:
fid = fopen([outPath,'/importBackgroundStimParams.hoc'],'w');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref BG_Stim_Dend[%d][%d], BG_Stim_Apic[%d][%d], BG_Stim_PV[%d][%d], BG_Stim_CB[%d][%d], BG_Stim_CR[%d][%d]\n',...
    size(bgStimInputDend,1),size(bgStimInputDend,2),size(bgStimInputApic,1),size(bgStimInputApic,2),size(bgStimInputPV,1),size(bgStimInputPV,2),...
    size(bgStimInputCB,1),size(bgStimInputCB,2),size(bgStimInputCR,1),size(bgStimInputCR,2));
fprintf(fid,'BG_dendSyn = %d\n',size(bgStimInputDend,2));
fprintf(fid,'BG_apicSyn = %d\n',size(bgStimInputApic,2));
fprintf(fid,'BG_PVSyn = %d\n',size(bgStimInputPV,2));
fprintf(fid,'BG_CBSyn = %d\n',size(bgStimInputCB,2));
fprintf(fid,'BG_CRSyn = %d\n',size(bgStimInputCR,2));

fprintf(fid,'\n\n// Import parameters: \n\n');
% Only for Pyramidals:
for c=1:size(bgStimInputDend,1)
    for i=1:size(bgStimInputDend,2)
        fprintf(fid,'BG_Stim_Dend[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputDend{c,i}));
        for j=1:length(bgStimInputDend{c,i})
            fprintf(fid,'BG_Stim_Dend[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputDend{c,i}(j)));
        end
    end
end


for c=1:size(bgStimInputApic,1)
    for i=1:size(bgStimInputApic,2)
        fprintf(fid,'BG_Stim_Apic[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputApic{c,i}));
        for j=1:length(bgStimInputApic{c,i})
            fprintf(fid,'BG_Stim_Apic[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputApic{c,i}(j)));
        end
    end
end


for c=1:size(bgStimInputPV,1)
    for i=1:size(bgStimInputPV,2)
        fprintf(fid,'BG_Stim_PV[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputPV{c,i}));
        for j=1:length(bgStimInputPV{c,i})
            fprintf(fid,'BG_Stim_PV[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputPV{c,i}(j)));
        end
    end
end

for c=1:size(bgStimInputCB,1)
    for i=1:size(bgStimInputCB,2)
        fprintf(fid,'BG_Stim_CB[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputCB{c,i}));
        for j=1:length(bgStimInputCB{c,i})
            fprintf(fid,'BG_Stim_CB[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputCB{c,i}(j)));
        end
    end
end
for c=1:size(bgStimInputCR,1)
    for i=1:size(bgStimInputCR,2)
        fprintf(fid,'BG_Stim_CR[%d][%d] = new Vector(%d)\n',c-1,i-1,length(bgStimInputCR{c,i}));
        for j=1:length(bgStimInputCR{c,i})
            fprintf(fid,'BG_Stim_CR[%d][%d].x[%d] = %d\n',c-1,i-1,j-1, abs(bgStimInputCR{c,i}(j)));
        end
    end
end
fclose(fid);
end