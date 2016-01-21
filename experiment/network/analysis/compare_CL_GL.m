pathprefix = 'C:/Users/user/Desktop/TEMP';

%%
run.nruns=1;
run.NC_str = 1;
run.nPC = 75;
run.ISBINARY = 1;
run.sn=12;
run.state=1;
%%

RUNS_str = {};
Sid=1;
fprintf('Loading runs...');
PCcells_str = {};
for stc=1:run.NC_str(Sid)
    stc
    for ru = 1:run.nruns
        for c=1:run.nPC
            if (run.ISBINARY)
                mycell = ncell(nrn_vread(sprintf('%s/STR_SN%d_ST%d_CL/%d_%d_%d.bin',pathprefix,run.sn,run.state,stc-1,c-1,ru-1),'n'),10);
            else
                mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,stc-1,c-1,ru-1)),10);
            end
            PCcells_str{c,ru} = mycell ;
        end
    end
    RUNS_str_CL{1,stc} = PCcells_str(:,:);
end
fprintf('DONE!\n');

%%

RUNS_str = {};
Sid=1;
fprintf('Loading runs...');
PCcells_str = {};
for stc=1:run.NC_str(Sid)
    stc
    for ru = 1:run.nruns
        for c=1:run.nPC
            if (run.ISBINARY)
                mycell = ncell(nrn_vread(sprintf('%s/STR_SN%d_ST%d_GL/%d_%d_%d.bin',pathprefix,run.sn,run.state,stc-1,c-1,ru-1),'n'),10);
            else
                mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,stc-1,c-1,ru-1)),10);
            end
            PCcells_str{c,ru} = mycell ;
        end
    end
    RUNS_str_GL{1,stc} = PCcells_str(:,:);
end
fprintf('DONE!\n');

%%
plot(RUNS_str_CL_S{1,stc}{1}.mv,'k');hold on;
plot(RUNS_str_CL{1,stc}{1}.mv,'b');hold on;
plot(RUNS_str_GL{1,stc}{1}.mv,'r')