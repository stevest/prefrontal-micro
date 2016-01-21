clc;
unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
tic;
t=1; % Initial without vclamps
VCLAMP = 0;
EXPERIMENT = 0;  % 0 = random, 1= structured gia to $#$#% NEURON
SIMPLIFIED = 1;
PARALLEL = 1;
if (SIMPLIFIED)
    neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
    '%smechanism_simple/x86_64/special -mpi -nobanner '...
    '-c "PARALLEL=%d" '...
    '-c "SIMPLIFIED=%d" '...
    '-c "CLUSTER_ID=%d" '...
    '-c "EXPERIMENT=%d" '...
    '-c "AA=%d" '...
    '-c "VCLAMP=%f" '...
    'final.hoc'], ...
    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,aa, VCLAMP);
else
    neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
    '%smechanism_complex/x86_64/special -mpi -nobanner '...
    '-c "PARALLEL=%d" '...
    '-c "SIMPLIFIED=%d" '...
    '-c "CLUSTER_ID=%d" '...
    '-c "EXPERIMENT=%d" '...
    '-c "AA=%d" '...
    '-c "VCLAMP=%f" '...
    'final.hoc'], ...
    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,aa, VCLAMP);
end
[nrnStatus,nrnCmdOut{t}] = unix(neuroncmd);
runtime = toc

% %
% clc
% nrnCmdOut{t}

%% Plot
plot(RUNS_str{1,t}{13,ru}.mv);
%% Plot
% close all;

ru = 1 ;% what run of this cluster

figure();
for c=1:nPC
    c
    comprobj(c) = RUNS_str{1,t}{c,ru};
    FFS(c) = RUNS_str{1,t}{c,ru}.spike_count(1,1500) / 3.5;
    FFPA(c) = RUNS_str{1,t}{c,ru}.spike_count(1500,tstop) / 3.5;
end
plot(FFS);hold on;

tmp=find(StimVect_str);
for c = 1:length(tmp)
plot(tmp(c),RUNS_str{1,t}{tmp(c),ru}.spike_count(1500,tstop) / 3.5,'+r');hold on;
end

