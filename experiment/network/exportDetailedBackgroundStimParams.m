function exportDetailedBackgroundStimParams(bgStimInputDend,bgStimInputApicpr,bgStimInputApic,bgStimInputSomaPV,bgStimInputSomaCB,bgStimInputSomaCR,outPath)

fid = fopen([outPath,'/importBackgroundStimParams.hoc'],'w');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref BG_Stim_basal[%d][%d][%d], BG_Stim_Apicpr[%d][%d][%d], BG_Stim_Apic[%d][%d][%d], BG_Stim_SomaPV[%d][%d][%d], BG_Stim_SomaCB[%d][%d][%d], BG_Stim_SomaCR[%d][%d][%d]\n',...
    size(bgStimInputDend,1),size(bgStimInputDend,2),size(bgStimInputDend,3),...
    size(bgStimInputApicpr,1),size(bgStimInputApicpr,2),size(bgStimInputApicpr,3),...
    size(bgStimInputApic,1), size(bgStimInputApic,2),size(bgStimInputApic,3),...
    size(bgStimInputSomaPV,1),size(bgStimInputSomaPV,2),size(bgStimInputSomaPV,3),...
    size(bgStimInputSomaCB,1),size(bgStimInputSomaCB,2),size(bgStimInputSomaCB,3),...
    size(bgStimInputSomaCR,1),size(bgStimInputSomaCR,2),size(bgStimInputSomaCR,3));

fprintf(fid,'BG_dendSyn = %d\n',size(bgStimInputDend,3));
fprintf(fid,'BG_apicSyn = %d\n',size(bgStimInputApic,3));
fprintf(fid,'BG_apicprSyn = %d\n',size(bgStimInputApicpr,3));
fprintf(fid,'BG_PVSyn = %d\n',size(bgStimInputSomaPV,3));
fprintf(fid,'BG_CBSyn = %d\n',size(bgStimInputSomaCB,3));
fprintf(fid,'BG_CRSyn = %d\n',size(bgStimInputSomaCR,3));



fprintf(fid,'\n\n// Import parameters: \n\n');
for many_times=1:size(bgStimInputDend,1)
for c=1:size(bgStimInputDend,2)
    for i=1:size(bgStimInputDend,3)
        fprintf(fid,'BG_Stim_basal[%d][%d][%d] = new Vector(%d)\n',many_times-1, c-1,i-1,length(bgStimInputDend{many_times, c,i}));
        for j=1:length(bgStimInputDend{many_times, c,i})
            fprintf(fid,'BG_Stim_basal[%d][%d][%d].x[%d] = %d\n',many_times-1, c-1,i-1,j-1, abs(bgStimInputDend{many_times,c,i}(j)));
        end
    end
end
end

for many_times=1:size(bgStimInputApicpr,1)
for c=1:size(bgStimInputApicpr,2)
    for i=1:size(bgStimInputApicpr,3)
        fprintf(fid,'BG_Stim_Apicpr[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(bgStimInputApicpr{many_times, c,i}));
        for j=1:length(bgStimInputApicpr{many_times, c,i})
            fprintf(fid,'BG_Stim_Apicpr[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(bgStimInputApicpr{many_times, c,i}(j)));
        end
    end
end
end

for many_times=1:size(bgStimInputApic,1)
for c=1:size(bgStimInputApic,2)
    for i=1:size(bgStimInputApic,3)
        fprintf(fid,'BG_Stim_Apic[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(bgStimInputApic{many_times, c,i}));
        for j=1:length(bgStimInputApic{many_times, c,i})
            fprintf(fid,'BG_Stim_Apic[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(bgStimInputApic{many_times, c,i}(j)));
        end
    end
end
end

for many_times=1:size(bgStimInputSomaPV,1)
for c=1:size(bgStimInputSomaPV,2)
    for i=1:size(bgStimInputSomaPV,3)
        fprintf(fid,'BG_Stim_SomaPV[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(bgStimInputSomaPV{many_times, c,i}));
        for j=1:length(bgStimInputSomaPV{many_times, c,i})
            fprintf(fid,'BG_Stim_SomaPV[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(bgStimInputSomaPV{many_times, c,i}(j)));
        end
    end
end
end

for many_times=1:size(bgStimInputSomaCB,1)
for c=1:size(bgStimInputSomaCB,2)
    for i=1:size(bgStimInputSomaCB,3)
        fprintf(fid,'BG_Stim_SomaCB[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(bgStimInputSomaCB{many_times, c,i}));
        for j=1:length(bgStimInputSomaCB{many_times, c,i})
            fprintf(fid,'BG_Stim_SomaCB[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(bgStimInputSomaCB{many_times, c,i}(j)));
        end
    end
end
end

for many_times=1:size(bgStimInputSomaCR,1)
for c=1:size(bgStimInputSomaCR,2)
    for i=1:size(bgStimInputSomaCR,3)
        fprintf(fid,'BG_Stim_SomaCR[%d][%d][%d] = new Vector(%d)\n',many_times-1,c-1,i-1,length(bgStimInputSomaCR{many_times, c,i}));
        for j=1:length(bgStimInputSomaCR{many_times, c,i})
            fprintf(fid,'BG_Stim_SomaCR[%d][%d][%d].x[%d] = %d\n',many_times-1,c-1,i-1,j-1, abs(bgStimInputSomaCR{many_times, c,i}(j)));
        end
    end
end
end

fclose(fid);
end