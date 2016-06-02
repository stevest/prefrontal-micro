
%% THis is the nrun_example.m from GIT. so put them tidely in the repo please...
% Load a BIG network and generate background activity:
cd('C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network')
load('states_07_G.mat')
pathprefix = 'Z:/data/GliaBackup/'
sn = 16
ID = 12;
SN = sn;
ST = 7;
%% Batch generate background activity!!! pou se xanw pou se briskw!
pathprefix='X:\Documents\Glia\'
for fores = 97:-1:50
tevents = cell(evlen,1);
% snormal = @(x,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,0,sigma);
kernel = snormal(-10:0.01:10,0,10,0.9);
sevents = cell(run.nruns,1);
rate = zeros(evlen,tstop);
% figure;hold on;
% we want the same average FF (1Hz) but 
events = 100;
for kk=1:evlen % input neurons BG
    kk/evlen
    tmp = zeros(1,tstop);
    rndshuffle = rand(1,events)*((tstop/events)*1); %the factor changes Vm variability!
    rndshuffle = rndshuffle - ((tstop/events)*1)/2;
    rndshuffle = round(linspace(10,tstop,events) + rndshuffle);
    rndshuffle(rndshuffle<=0) = [];
    tmp(rndshuffle) = 10;
    rate(kk,1:tstop) = conv(tmp(1,1:tstop),kernel,'same');
%     tmp = rate(kk,1:tstop);    
tevents{kk} = nonhomogenouspoissrnd(tstop,0.12,rate(kk,:));
% plot(tevents{kk},1,'*')
%     run.fastSpikesMatBackground = cell2mat(sevents(1,1));
%     run.exportBackgroundStimParams(pathprefix)
%     movefile(sprintf('%s/importBackgroundStimParams_run_0001.hoc',pathprefix),sprintf('%s/importBackgroundStimParams_run_%04d.hoc',pathprefix,kk));
end
maxLength = max(cellfun(@length,tevents));
sevents{1,1}=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))*NaN),tevents,'UniformOutput',false));
run.fastSpikesMatBackground = cell2mat(sevents(1,1));

run.exportBackgroundStimParams(pathprefix)
movefile(sprintf('%simportBackgroundStimParams_run_0001.hoc',pathprefix),sprintf('%simportBackgroundStimParams_run_0%03d.hoc',pathprefix,fores))
end
%% LOAD THE RESPECTIVE BIG NETWORK DATA

% override states:
states.gparams(7).nPC = nPC;
states.gparams(7).nPV = nPV;
states.gparams(7).nCB = 1;
states.gparams(7).nCR = 1;
states.state_str  = PC2PC_str(:,:,:);
states.state_rnd  = PC2PC_rnd(:,:,:);
%Fill the rest states if not filled:
states.state_str  = repmat(PC2PC_str(:,:,:),[1,1,10]);
states.state_rnd  = repmat(PC2PC_rnd(:,:,:),[1,1,10]);

run = nrun(ID,nPC,SN,ST,1,10000);
%     run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
run.pathToHere = 'C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network';
run.SIMPLIFIED = 1;
% Insert cynaptic clustering bias: 0 = clustered, 1 = random (no
% clustering)
run.CLUSTBIAS = 0; % keep synapses in original locations (no rand jittering)
run.FORCECLUSTERING = 0; % Recompute the clusters (overriding the precomputed states)
run.init(states);

tstop = run.tstop;
evlen = 20;
tevents = cell(evlen,1);
snormal = @(x,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,0,sigma);
kernel = snormal(-10:0.01:10,10,0.9);
sevents = cell(run.nruns,1);
rate = zeros(evlen,tstop);
% figure;hold on;
% we want the same average FF (1Hz) but 
events = 10;
for kk=1:evlen % input neurons BG
    tmp = zeros(1,tstop);
    rndshuffle = rand(1,events)*((tstop/events)*1); %the factor changes Vm variability!
    rndshuffle = rndshuffle - ((tstop/events)*1)/2;
    tmp(round(linspace(10,tstop,events) + rndshuffle)) = 10;
    rate(kk,1:tstop) = conv(tmp(1,1:tstop),kernel,'same');
    tmp = rate(kk,1:tstop);    
    % compute rate variability (approximate):
    nrate = tmp./(max(tmp)/0.12);
%     plot(histc(nrate,0:0.01:1))
    std(tmp)
    
tevents{kk} = nonhomogenouspoissrnd(tstop,0.12,rate(kk,:));
% plot(tevents{kk},1,'*')
%     run.fastSpikesMatBackground = cell2mat(sevents(1,1));
%     run.exportBackgroundStimParams(pathprefix)
%     movefile(sprintf('%s/importBackgroundStimParams_run_0001.hoc',pathprefix),sprintf('%s/importBackgroundStimParams_run_%04d.hoc',pathprefix,kk));
end
% create lognormal weights:
% lognormalW = lognpdf(0:0.01:2.3,0.49,5);
% figure;plot(0:0.01:2.3,lognormalW)
lognormalW = lognrnd(0.49,5,[1,evlen]);

figure;plot(sum(rate))
figure;plot(rate' * lognormalW')
%     maxLength = max(cellfun(@length,tevents));
%     sevents{1,1}=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))*NaN),tevents,'UniformOutput',false));
    
tmp = rate(kk,1:tstop);
Fs = tstop;
N = length(tmp);
xdft = fft(tmp);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(tmp):Fs/2;
figure;
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
set(gca,'xscale','log')

%% Previous RunExample
clear all;close all;clc;
cd('C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network')
% cd('E:\NEURON_RUNS')
load('states_07_G.mat')
% pathprefix = 'N:/NEURON_PROJECTS/NEW_RUNS/';
% pathprefix = 'H:/NEURON_RUNS/';
% pathprefix = 'I:/data/demory_backup/NEURON_PROJECTS/NEW_RUNS/';
% pathprefix = 'C:/Users/user/Desktop/TEMP/TEMP/';
% pathprefix = 'E:/NEURON_RUNS/';
pathprefix = 'Z:/data/GliaBackup/'
% pathprefix = 'Z:/data/demory_backup/NEURON_PROJECTS/NEW_RUNS/';
% pathprefix = 'Z:/data/demory_backup/NEURON_PROJECTS/NEW_RUNS/TEMP/'

% states.PC2PC_rnd = PC2PC_rnd;
% states.PC2PC_str = PC2PC_str;

% for currentstate = 1:11
%     run = nrun(12,75,1,currentstate,10,5000);
%     run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
%     run.init(states);
%     states.Labels_str(:,currentstate) = run.labels_str;
%     states.Labels_rnd(:,currentstate) = run.labels_rnd;
% end

    
%%
% clustbiasRange = 0:0.1:1;
for sn = 16
    clearvars -except states sn tmpSPK pathprefix clustbiasRange
    % obj = nrun(experimentid,npyrs,serial,state,exprun,tstop)
    ID = 12;
    SN = sn;
    ST = 7;
    % Manipulate the STR case to get similar GCC and No of reciprocals as in the rnd case:
    p=0.5;
    ch = states.state_str(1:75,1:75,ST);
    for k=1:75
        for l=1:75
            if (ch(k,l) && (rand<p) )
                tmpk=1;tmpl=1;
                while (tmpk == tmpl) || (ch(tmpk,tmpl)) % gia to deftero conditional den eimai sigouros.
                    [tmpk,tmpl] = ind2sub(size(ch),floor(1+rand*75*75));
                end
                ch(tmpk,tmpl) = 1;
                ch(k,l) = 0;
            end
        end
    end
    tmpstate = states.state_str(1:75,1:75,ST);
    [GCC_str,LCC_str,~] = clust_coeff( tmpstate );
    fprintf('CC:%f Total:%f Recurent:%2.3f%%\n', GCC_str, sum(tmpstate(:)), 100*sum(sum(tmpstate & tmpstate'))/sum(tmpstate(:)))
    [mcccs_ch,lcc_ch,~] = clust_coeff( ch );
    fprintf('CC:%f Total:%f Recurent:%2.3f%%\n', mcccs_ch, sum(ch(:)), 100*sum(sum(ch & ch'))/sum(ch(:)))
    % Overwrite changes in states variable!
    states.state_str(1:75,1:75,ST) = ch;
    
    

    run = nrun(ID,1000,SN,ST,1,10000);
%     run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
    run.pathToHere = 'C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network';
    run.SIMPLIFIED = 1;
        % Insert cynaptic clustering bias: 0 = clustered, 1 = random (no
    % clustering)
    run.CLUSTBIAS = 0; % keep synapses in original locations (no rand jittering)
    run.FORCECLUSTERING = 0; % Recompute the clusters (overriding the precomputed states)
    run.init(states);
        
    
%     sum(sum(run.stateSTR(1:run.nPC,1:run.nPC)))            / (run.nPC*run.nPC)
%     mkdir(run.path);
%     save(sprintf('%s/EXP_ID%d_SN%d_ST%d.mat',run.path,run.id,run.sn,run.state),'run');

    
    
%     for i=1:run.NC_str
%         
%         figure(1);
%         i
%         ci = [];
%         ci = find(run.labels_str==i);
%         tmpmat = run.stateSTR(ci,ci);
%         s = ismember(ci,run.StimMat_str(:,i));
%         %     tmpmat(s,s) = 2;
%         sc(tmpmat);hold on;
%         scatter(ones(1,7),find(s),70,'r','o','fil')
%         connPerClust(i) = sum(sum(run.stateSTR(ci,ci)))
%         %     pause;
%         %     cla;
%         close(1)
%     end
%     
%     [~,mostCellsCluster] =max(run.cellsPerCluster_str)
%     [~,maxCellsCluster] = max(connPerClust);
%     thisCluster = maxCellsCluster;
    
    % run(clusterid,EXPERIMENT,SIMPLIFIED,PARALLEL)
    
    % Export Parameters to NEURON:
    run.exportParams2Neuron(pathprefix);
%     run.exportStimulationParameters(pathprefix)
%     run.exportNetworkParameters(pathprefix)
%     run.exportNetworkStimulationHeader(pathprefix)
%     run.exportNetworkStimulation(pathprefix)
%     run.exportBackgroundStimParamsHeader(pathprefix)
%     run.exportBackgroundStimParams(pathprefix)

tstop = run.tstop;
evlen = size(run.fastSpikesMatBackground,1)/run.nruns;
events = cell(run.nruns,1);
for kk=1:1%run.nruns
    rate = zeros(1,tstop);
    rate(round([10:1000:tstop] + rand(1,10)*700)) = 10;
    % rate(find(poissrnd(0.002,1,tstop))) = 10;
    % create skewed Normal Distribution as BA:
    snormal = @(x,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,0,sigma);
    % kernel = fspecial('gaussian',[1 1000], 250/3);
    kernel = snormal(-10:0.01:10,10,0.9);
%     figure;plot(kernel)
    rate = conv(rate,kernel,'same');
    % plot(rate);

    % events = zeros(8752,tstop);
    tevents = cell(evlen,1);
%     figure;hold on;
    for k=1:evlen
        k
        tevents{k} = nonhomogenouspoissrnd(tstop,0.006,rate);
    %     events(k,tevents) = 1;
    %     scatter(events(k,:),ones(1,length(events))*k,'+');
    end
    % plot(sum(events))

    maxLength = max(cellfun(@length,tevents));
    events{1,1}=cell2mat(cellfun(@(x)cat(2,x,zeros(1,maxLength-length(x))*NaN),tevents,'UniformOutput',false));
    
    run.fastSpikesMatBackground = cell2mat(events(1,1));
    run.exportBackgroundStimParams(pathprefix)
    movefile(sprintf('%s/importBackgroundStimParams_run_0001.hoc',pathprefix),sprintf('%s/importBackgroundStimParams_run_%04d.hoc',pathprefix,kk));
end

run.fastSpikesMatBackground = cell2mat(events);
run.exportBackgroundStimParamsHeader(pathprefix)
run.exportBackgroundStimParams(pathprefix)
% Visualize background activity pontiako implementation:
viz = zeros(round(size(run.fastSpikesMatBackground,1)/100),100);
for k=1:size(run.fastSpikesMatBackground,1)/100
    viz(k,:) = histcounts(run.fastSpikesMatBackground(k,:),0:100:run.tstop);
end
figure;hold on;
imagesc(viz);

    mkdir(pathprefix,run.path);
    run.fastSpikesMatStimulation = [];
    run.fastSpikesMatBackground = [];
    run.SpikesStimDend = {};
    run.SpikesStimApic = {};
    run.SpikesBgBasal = {};
    run.SpikesBgProximal = {};
    run.SpikesBgApical = {};
    run.SpikesBgPV = {};
    run.SpikesBgCB = {};
    run.SpikesBgCR = {};
    
    save(sprintf('%s%s/EXP_ID%d_SN%d_ST%d.mat',pathprefix,run.path,run.id,run.sn,run.state),'run','-v7.3');

    % Push data to cluster:
%     run.push();
    
    whatnet = 1; % 0 = random, 1= structured gia to $#$#% NEURON
    if whatnet
        nclu = run.NC_str;
    else
        nclu = run.NC_rnd;
    end
                
    for ccl = 1:nclu
        run.run(ccl,whatnet,1,1);
    end
    
            

end
% pull results from cluster:
% run.pull();

%% Load STR
% close all; clear all; clc;
% load(sprintf('STR_%d_%d.mat',12,8)) % 8,10
% load(sprintf('%sexperiment_%d/EXP_ID%d_SN%d_ST%d.mat',pathprefix,12,12,15,1)) % 8,10
% mat = matfile(sprintf('%sexperiment_%d/EXP_ID%d_SN%d_ST%d.mat',pathprefix,12,12,15,1)) % 8,10

RUNS_str = cell(1,run.NC_str);
Sid=1;
fprintf('Loading runs...');
PCcells_str = cell(run.nPC,run.nruns);
ridx = zeros(run.NC_str(Sid),run.nPC,run.nruns);
for stc=1%:run.NC_str(Sid)
    RUNS_str{1,stc} = cell(run.nPC,run.nruns);
    for ru = 1:run.nruns
        ru
        tic;
        S = cell(1,run.nPC);
        parfor pc=1:run.nPC
            if (run.ISBINARY)
                try 
                    PCcells_str{pc,ru} = ncell(nrn_vread(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1),'n'),10);
                catch err
                    warning('on','verbose')
                    warning(err.message)
                    ridx(stc,pc,ru) = 1; % true, if file not found
                    continue;
                end
            else
                % Following line (load) not working with parfor:
%                 mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,stc-1,pc-1,ru-1)),10);
            end
            %             mycell = ncell(load(sprintf('%sexperiment_10_10Hz_Stim/%s/%d_%d_%d.txt',mypath,'STR',t,c-1,ru-1)),10);
            PCcells_str{pc,ru}.clusterID = run.labels_str(pc,Sid);
            %             mycell.position = PCsomata(c,1:3);
            %             PCcells_str{c,ru}=mycell.hasPersistent(run.stimend-1,25,run.tstop-(run.stimend-1)); % paper?
            [S{pc},~,~] = findUPstates(PCcells_str{pc,ru}.mv(run.stimend*run.dt:run.dt:end),4, 10, -66, 3000 );
%             if ~isempty(S)
%                 plot(mycell.mv)
%                 pause();
%                 cla;
%             end
            % For MEMORY concerns reduce object's size:
%             mycell.mv = single(mycell.mv(1:mycell.dt:end));
%             mycell.dt=1;
            PCcells_str{pc,ru}.mv = [];
            PCcells_str{pc,ru}.spikes = single(PCcells_str{pc,ru}.spikes);
            
            if ~isempty(S{pc});
                PCcells_str{pc,ru}.persistentActivity = 1;
            else
                PCcells_str{pc,ru}.persistentActivity = 0;
            end
            PCcells_str{pc,ru} = PCcells_str{pc,ru} ;
            %             leastSynapses(c) = PCcells_str{c,ru}.nspikes;
        end
        runtime = toc - tic
%                 for c=run.nPC+1:run.nPC+run.nPV
%                     mycell = ncell(nrn_vread(sprintf('%s%s/STR_%d/%d_%d_%d.bin',pathprefix,run.path,run.state,clu-1,c-1,ru-1),'n'),10);
%                     PVcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
%                 end
        %         for c=run.nPV+1:run.nPV+run.nCB
        %             mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,t-1,c-1,ru-1)),10);
        %             CBcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        %         end
        %         for c=run.nCB+1:run.nCB+run.nCR
        %             mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,t-1,c-1,ru-1)),10);
        %             CRcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        %         end
    end
    RUNS_str{1,stc} = PCcells_str(:,:);
    % rescale the cell to discard empty runs (with no file loaded):
%     for k = 1:size(RUNS_str{1,stc},2)
%         ridx(stc,k) = isempty(RUNS_str{1,stc}{1,k}.spikes);
%     end
%     RUNS_str{1,stc}(:,find(ridx(stc,:))) = [];
    RUNS_str{1,stc}(:,find(any(ridx(stc,:,1:size(PCcells_str(:,:),2)),2))) = [];
end

% run.nruns = min(sum(~ridx,2));
% run.nruns = size(RUNS_str{1,stc},2)

% Save memory!
clear PCcells_str;

fprintf('DONE!\n');
%% Load RND
% close all; clear all; clc;
% load(sprintf('STR_%d_%d.mat',12,2))
% load(sprintf('%sexperiment_%d/EXP_ID%d_SN%d_ST%d.mat',pathprefix,12,12,2,1)) % 8,10

RUNS_rnd = cell(1,run.NC_rnd);
Sid=1;
fprintf('Loading runs...');
PCcells_rnd = {};
ridx = zeros(run.NC_rnd(Sid),run.nPC,run.nruns);
for stc=1%:run.NC_str(Sid)
    RUNS_rnd{1,stc} = cell(run.nPC,run.nruns);
    for ru = 1:run.nruns
        ru
        for c=1:run.nPC
            if (run.ISBINARY)
                try 
                    mycell = ncell(nrn_vread(sprintf('%s%s/RND_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,c-1,ru-1),'n'),10);
                catch err
                    warning('on','verbose')
                    warning(err.message)
                    ridx(stc,c,ru) = 1; % true, if file not found
                    continue;
                end
            else
                mycell = ncell(load(sprintf('%s/RND_%d/%d_%d_%d.txt',run.path,run.state,stc-1,c-1,ru-1)),10);
            end

            mycell.clusterID = run.labels_rnd(c,Sid);

            [S,~,~] = findUPstates(mycell.mv(run.stimend*run.dt:run.dt:end),4, 10, -66, 3000 );

            mycell.mv = [];
            mycell.spikes = single(mycell.spikes);
            
            if ~isempty(S);
                mycell.persistentActivity = 1;
            else
                mycell.persistentActivity = 0;
            end
            PCcells_rnd{c,ru} = mycell ;
            %             leastSynapses(c) = PCcells_str{c,ru}.nspikes;
        end
    end
    RUNS_rnd{1,stc} = PCcells_rnd(:,:);
    RUNS_rnd{1,stc}(:,find(any(ridx(stc,:,1:size(PCcells_rnd(:,:),2)),2))) = [];
end

% Save memory!
clear PCcells_rnd;
fprintf('DONE!\n');

%%
figure();
for i=1:run.nPV
    plot(PVcells_str{75+i,1}.mv);hold on;
end
% plot(CBcells_str{89,1}.mv)
% plot(CRcells_str{95,1}.mv)
%% ammount of PV connections:
PVconns = run.stateSTR(run.nPC+1:run.nPC+run.nPV,1:run.nPC);
sum(sum(PVconns)) / (run.nPV*run.nPC)

%% Plot-a-cluster STR
close all;

clu= 1;% what cluster is stimulated
ru = 2 ;% what run of this cluster
Sid = 1;
CellsPerClust=zeros(run.NC_str(Sid),1);
PAPerClust=zeros(run.NC_str(Sid),1);


for cl=1:run.NC_str(Sid)
    figure(); hold on;
    set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
    for c=1:run.nPC
        if(RUNS_str{1,clu}{c,ru}.clusterID == cl)
            %                         figure;
            plot( RUNS_str{1,clu}{c,ru}.mv,'linewidth',1,'color',rand(1,3));
            CellsPerClust(cl) = CellsPerClust(cl)+1;
%             RUNS_str{1,clu}{c,ru} = RUNS_str{1,clu}{c,ru}.hasPersistent(run.stimend-1,60,run.tstop-(run.stimend-1)); % paper?
%             [S,D,UPs]=findUPstates(RUNS{1,clu}{c,ru}.mv(run.stimend:run.dt:end), 2, 3, -66, 3000 ) ;
            if(RUNS_str{1,clu}{c,ru}.persistentActivity)
                PAPerClust(cl) = PAPerClust(cl) +1;
            end
        end
        %         pause
    end
    axis([0,length(RUNS_str{1,clu}{c,ru}.mv),-70,70 ]);hold on;
end

% bar graph
figure();
bar(CellsPerClust,'b');hold on;
bar(PAPerClust,'r');
%% Plot-a-cluster RANDOM!
close all;

clu= 1;% what cluster is stimulated
ru =10 ;% what run of this cluster
Sid = 1;
CellsPerClust=zeros(run.NC_rnd(Sid),1);
PAPerClust=zeros(run.NC_rnd(Sid),1);


for cl=1:run.NC_rnd(Sid)
    figure(); hold on;
    set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
    for c=1:run.nPC
        if(RUNS_rnd{1,clu}{c,ru}.clusterID == cl)
            %                         figure;
            plot( RUNS_rnd{1,clu}{c,ru}.mv,'linewidth',1,'color',rand(1,3));
            CellsPerClust(cl) = CellsPerClust(cl)+1;
%             RUNS_str{1,clu}{c,ru} = RUNS_str{1,clu}{c,ru}.hasPersistent(run.stimend-1,60,run.tstop-(run.stimend-1)); % paper?
%             [S,D,UPs]=findUPstates(RUNS{1,clu}{c,ru}.mv(run.stimend:run.dt:end), 2, 3, -66, 3000 ) ;
            if(RUNS_rnd{1,clu}{c,ru}.persistentActivity)
                PAPerClust(cl) = PAPerClust(cl) +1;
            end
        end
        %         pause
    end
    axis([0,length(RUNS_rnd{1,clu}{c,ru}.mv),-70,70 ]);hold on;
end

% bar graph
figure();
bar(CellsPerClust,'b');hold on;
bar(PAPerClust,'r');

figure;
figure;
for i=76:88
    plot(PV_rnd{1,1}{i}.mv);hold on;
end

%% Calculate the stationary distribution of the network:
stc = 1;

Q = 2:10;% Q = 100ms integration time of the cell (?)
for ql = 1:length(Q)
    if Q(ql) >= run.stimend
        error('Q window length is too big!');
    end
    %unique states are N = 2^n.
    N = 2^run.stimCellsPerCluster;
    sProbs{ql} = zeros(N,length(run.stimend + Q(ql) : run.tstop));
    
    % get spike trains from the cells of interest:
    spktrain = zeros(run.stimCellsPerCluster,run.tstop,run.nruns);
    for ru = 1:run.nruns
        for c = 1:run.stimCellsPerCluster
            spktrain(c,round(RUNS_str{1,1}{run.StimMat_str(c,stc),ru}.spikes'),ru) = 1;
        end
    end
    
    starting = 2 + Q(ql) ;
    % starting = run.stimend + Q;
    % for t = run.stimend + Q :Q: run.tstop
    for t = starting:Q(ql): run.tstop
        
        % extract the states existing inside window in time t
        sstates = spktrain(:,(t-Q(ql))+1:t,:);
        %     sstates = zeros(run.stimCellsPerCluster,Q,run.nruns);
        %     for ru = 1:run.nruns
        %         for c = 1:run.stimCellsPerCluster
        %             spikePool = round(RUNS_str{1,stc}{run.StimMat_str(c,stc),ru}.spikes');
        %             sstates(c, unique(spikePool((spikePool <= t) & (spikePool > (t-Q)))) - (t-Q), ru) = 1 ; % unique in case we have miltiple spikes in a milisecond... kanonika den 8a eprepe na xreiazetai..
        %         end
        %     end
        sstates = reshape(any(sstates,2),size(any(sstates,2),1),[]);
        
        %     % Calculate probability for each state between runs:
        %     for l=1:size(sstates,2)
        %         histo(bi2de(sstates(:,l)')+1 , (t-(run.stimend + Q)) + 1 ) = histo(bi2de(sstates(:,l)')+1 , (t-(run.stimend + Q)) + 1 ) + 1;
        %     end
        [~,ia,ic] = unique(sstates','rows');
        for l=1:length(ia)
            sProbs{ql}(bi2de(sstates(:,ia(l))')+1 , (t-starting) + 1 ) = sum(ic == l) / run.nruns; % na dw mipws exw la8os edw..
        end
        
        fprintf('@ t=%d to result paei sto %d\n',t,(t-starting) + 1)
    end % for all the NEURON runs
end

% figure;sc(sProbs,jet);
figure;
bwin = 50;
for l=10
%     plot(smoothts(sProbs(l,:),'g',1000),'color',rand(1,3));hold on;
    csum = cumsum(sProbs(l,:));
    vara = bwin:bwin:size(sProbs,2);
    varb = 1:bwin:size(sProbs,2)-bwin;
    vari = csum(vara) - csum(varb);
    plot(smoothts(vari),'color',rand(1,3));hold on;
end


mcmcgr(sProbs{ql},ng)

figure	;
bwin = 50	;
l=10 ;
cm = jet(length(Q));
for ql=1:length(Q)
%     plot(smoothts(sProbs(l,:),'g',1000),'color',rand(1,3));hold on;
    csum = cumsum(sProbs{ql}(l,:));
    vara = bwin:bwin:size(sProbs{ql},2);
    varb = 1:bwin:size(sProbs{ql},2)-bwin;
    vari = csum(vara) - csum(varb);
    plot(smoothts(vari),'color',cm(ql,:));hold on;
end


% %% Calculate the stationary distribution of the network:
% stc = 1;
% for ru = 1:run.nruns
%     % create spiketrain matrix:
%     spkMat{1,ru} = zeros(run.stimCellsPerCluster,run.tstop);
%     for c = 1:run.stimCellsPerCluster
%         spkMat{1,ru}(c,round(RUNS_str{1,stc}{run.StimMat_str(c,stc),ru}.spikes')) = 1 ;
%     end
%     
%     Q = 5; % Q = 100ms integration time of the cell (?)
%     t = run.stimend + Q;
%     if Q >= run.stimend
%         error('Q window length is too big!');
%     end
%     sstates{1,ru} = [];
%     while run.tstop - t
%         sstates{1,ru} = [sstates{1,ru} any(spkMat{1,ru}(:,t-Q:t),2)];
%         t = t+1;
%     end
%     %unique states are N = 2^n. 
%     N = 2^run.stimCellsPerCluster;
%     histo{1,ru} = zeros(1,N);
%     for l=1:size(sstates{1,ru},2)
%         histo{1,ru}(1,bi2de(sstates{1,ru}(:,l)')+1) = histo{1,ru}(1,bi2de(sstates{1,ru}(:,l)')+1) + 1;
%     end
%     % compute relative ferquency in time (? how i do such a thing?):
%     % for arbitrary network state:
%     for Nest = 1:N;
%         for t = 2:size(sstates{1,ru},2)%run.stimend + Q +1 : run.tstop
%             tdif = t - 1;
%             P{1,ru}(Nest,t) = sum(ismember(sstates{1,ru}(:,1:t)',de2bi(Nest-1,run.stimCellsPerCluster),'rows')) / (tdif * Nest);
%         end
%     end
% %     [~,ia,ic] = unique(sstates','rows');
% %     sProbs = zeros(1,length(ia));
% %     for l = 1:length(ia)
% %         sProbs(l) = sum(ic == l);
% %     end
%     
% end % for all the NEURON runs

%% MCMCGR - Gelman-Rubin R statistic for convergence
% save(sprintf('%s%s/RUNS_str_ID%d_SN%d_ST%d_STC_%d.mat',pathprefix,run.path,run.id,run.sn,run.state,stc),'RUNS_str','-v7.3');
% save(sprintf('%s%s/RUNS_rnd_ID%d_SN%d_ST%d_STC_%d.mat',pathprefix,run.path,run.id,run.sn,run.state,stc),'RUNS_rnd','-v7.3');
% load(sprintf('%s%s/RUNS_str_ID%d_SN%d_ST%d_STC_%d_withBG.mat',pathprefix,run.path,run.id,run.sn,run.state,stc));
% load(sprintf('%s%s/RUNS_rnd_ID%d_SN%d_ST%d_STC_%d.mat',pathprefix,run.path,run.id,run.sn,run.state,stc));
RUNS = RUNS_rnd;
NC = run.NC_rnd;
% labels = run.labels_str;
labels = zeros(700,1);
for k = 1:run.NC_rnd
    labels(run.StimMat_rnd(:,k)) = k;
end
% save('RUNS_MSRD_STR','RUNS','-v7.3');
% get vector with valid runs:

% % If forgot to save failedToLoad, recreate it!
% failedToLoad = zeros(run.nPC,100,7);
% for stc=1:7
%     for r = 1:100
%         for c = 1:run.nPC
%             failedToLoad(c,r,stc) = isempty(RUNS{1,stc}{c,r});
%         end
%     end
% end
rurange = find(~any( squeeze(any(failedToLoad,1)),2 ));
for stc=1:7
        RUNS{1,stc} = RUNS{1,stc}(:,rurange);
end
run.nruns = length(rurange);

% Get probability for each state in a run:


fr = [];
fridx = [];
pa = [];
stc = 1
for ru = [1:run.nruns]
    for c = 1:run.nPC
        if isempty(RUNS{1,stc}{c,ru}.freq)
            fr(c,ru) = 0;
%             pa(c,ru) = 0;
            continue;
        end
        % manual frequency:
        fr(c,ru) = length(RUNS{1,stc}{c,ru}.spikes > 1500) ;
%         fr(c,ru) = RUNS{1,stc}{c,ru}.freq;
%         pa(c,ru) = RUNS{1,stc}{c,ru}.persistentActivity;
    end
end
cm = jet(ceil(max(fr(:))));
cm(1,:) = [0,0,0];
figure();imagesc(fr);

% figure();imagesc(fr(labels==1,:));
figure();imagesc(fr);
% [~,IDX] = sort(mean(fr(labels==1,:),2),'descend')
[~,IDX] = sort(mean(fr,2),'descend')
% tmp = find(labels==1);
tmp = 1:run.nPC;
cellpool = {};
% cellpool{1,stc} = tmp(IDX(1:100))
cellpool{1,stc} = tmp(IDX)
figure;title('Working only with cells:');
imagesc(fr(cellpool{1,stc},:));

%For all network!:
[~,IDX] = sort(mean(fr,2),'descend')
tmp = 1:run.nPC;
cellpool = {};
for kk = 2:2:length(tmp)
    cellpool = [cellpool,tmp(IDX(1:kk))];
end
% For stimulated cluster only!
figure();imagesc(fr(labels==1,:));
[~,IDX] = sort(mean(fr(labels==1,:),2),'descend')
tmp = find(labels==1);
cellpool = {};
cellpool{1,stc} = tmp(IDX(1:100))
figure;title('Working only with cells:');
imagesc(fr(cellpool{1,stc},:));



%For surf graph:
[~,IDX] = sort(mean(fr(labels==1,:),2),'descend')
tmp = find(labels==1);
cellpool = {};
% prepei na elegxw oti:
% a) den exw 0 activity sto mean twn Q tou kane batch
% b) den exw 0 activity sto parapanw, gia ka8e chain
for kk = 2:2:length(tmp)
    tmp2 = tmp(IDX(1:kk));
    tmp3 = all(any( wspktrain(tmp2,:,:), 3),2);
    cellpool = [cellpool,tmp2(tmp3)];
end
% Plot surf graph:
surfmat = [];
for kk = 1:size(cellpool,2)
    surfmat = [surfmat ; ALL_hat_R{1,1,kk}]
end
figure;surf(surfmat)


colormap(cm)
% figure();imagesc(pa);
[~, fridx ] = sort(fr(:,1));
%FF over time:
steps = run.dt;
bin = 100;
start = run.stimstart;
bins = start:bin:stop;
freqBinsCell = cell(1,50);
clstidx = 1:700;%run.StimMat_str(:);
for ru = rurange
    freqBins = zeros(run.nPC,length(bins));
    for stc = 1%:run.NC_str
        for k = 1:length(clstidx)
            spikes = RUNS{1,stc}{clstidx(k),ru}.spikes;
            freqhisto = (histc(spikes, bins) / (bin/1000));
            if ~isempty(freqhisto)
                freqBins(k,:) = freqhisto ;
            end
        end
    end
    freqBinsCell{ru} = freqBins;
end
logicFreq = sum(reshape(cell2mat(cellfun(@(x) logical(x), freqBinsCell,'uniformoutput',false)),run.nPC,length(bins),[]),3);
figure;imagesc(logicFreq(:,11:end))
[varCells,lessVarCells] = sort(var(logicFreq(:,11:end),0,2)) ;
zeroRuns = any(logicFreq,2);
figure;imagesc(logicFreq(zeroRuns,:))
% zeroRuns = lessVarCells( (varCells<1) & (varCells~=0) )
[maxCells,lessMaxCells] = sort(max(logicFreq(:,11:end),[],2)) ;

zeroRuns = lessMaxCells( (maxCells>40) )
cellpool{1,stc} = (zeroRuns);
figure;imagesc(freqBinsCell{2}(cellpool{1,stc},:))

cells = find(all(fr'));
% RUNS{1,stc}(:,find(~all(pa))) = [];
run.nruns = size(RUNS{1,stc},2)
run.tstop = 20000 ;
RUNS{1,stc}(:,34) = []
RUNS{1,stc}(:,15) = []


% get spike trains from the cells of interest:
% sp = floor(run.nPC / NC);
% rp = randperm(run.nPC);
% rp = reshape(rp(1:sp*NC),NC,[]);
for stc = 1:NC
%     cellpool{1,stc} = find(labels == stc)';

% different handling:
    cellpool{1,stc} = run.StimMat_str(:,stc);
    
    
% Get random cells to check if convergence is specific to a cluster :
%      cellpool{1,stc} = rp(stc);
end
% Choose the runs that the cluster #stc was stimulated:
stc = 1;
cellpool{1,stc} = find(all(fr'));


%%
% ALL THE CODE BELOW follow the NOTATION as in:
% General methods for monitoring convergence of iterative simulations,
% Brooks and Gelman, 1998 !
% In parentheses the variables are notated EXACTLY as in the paper!

m = run.nruns ; %No of chains (m)
n = tstop; %No of itterations (n)

% Qseq = [2,5,10,15,20,30,40,80,150,200,500];
% Qseq = [2,3,4,5,6,7];
Qseq = [30];
Nseq = 1;

ALL_hat_R = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_GR_B_N = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_GR_W = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_hat_V = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_bar_a = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;


for Qs = 1:length(Qseq) 
    for Ns = 1:length(Nseq) 
        Q = Qseq(Qs) ; % simple window (ms)
        Qr = floor(n / Q) ; % length of wspiketrain
        n_over_b = 20 ; % No of batches
        b = floor(Qr/n_over_b); % Batches' length (b)
        if b < 2
            error('Batches'' length must be greater than one sample!');
        end
        
%         dist_WV = zeros(size(cellpool,2),n_over_b) ;
        hat_R = zeros(size(cellpool,2),n_over_b) ;
        GR_W = cell(size(cellpool,2),n_over_b) ;
        GR_B_N = cell(size(cellpool,2),n_over_b) ;
        hat_V = cell(size(cellpool,2),n_over_b) ;
        bar_Xi = zeros(run.nPC,n_over_b,m) ;
        dd_X = zeros(run.nPC,n_over_b) ;
        bar_a = cell(size(cellpool,2),n_over_b);
        
        s = (((1:Qr)-1) * Q) + 1 ;
        e = ((1:Qr) * Q) ;
        ts = [];
        te = [];
        x_len = [];
        
        % Calculate spike trains ONLY for the above selected cells:
        wspktrain = zeros(run.nPC,Qr,m);
        for ru = 1:m
            spktrain = zeros(run.nPC,n);
            for c=1:run.nPC % must be row vector!!
                if ~isempty(RUNS{1,stc}{c,ru})
                    spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
                end
            end
            for k=1:Qr
                wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
            end
        end       
        
        % For each Q now, get probability (per run):
        %         figure;hold on;
        %         cm = hsv(m);
        % Find from all the posible states nPC^2, which are represented in my
        % sample:
        clusterStatesAll = cell(1,7);
        voteState = cell(1,7);
        for stc=1:NC
            % Calculate spike trains ONLY for the above selected cells:
            wspktrain = zeros(run.nPC,Qr,m);
            for ru = 1:m
                spktrain = zeros(run.nPC,n);
                for c=1:run.nPC % must be row vector!!
                    if ~isempty(RUNS{1,stc}{c,ru})
                        spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
                    end
                end
                for k=1:Qr
                    wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
                end
            end

            % Break to batches to avoid out of memory:
            chunks = 1:700*100:numel(wspktrain);
            clusterStates = {};
            for chunk = 1:length(chunks)-1
                [clusterStates{chunk}, ~, ~] = unique(reshape(wspktrain(chunks(chunk):chunks(chunk+1)-1),run.nPC,[])', 'rows');
            end
            clusterStatesUnion = [];
            clusterStatesUnion = cell2mat(cellfun(@(x)x',clusterStates,'uniformoutput',false));
            [clusterStatesAll{stc}, ~, ~] = unique(clusterStatesUnion', 'rows');
    %         % If no batches are used:
    %         [clusterStates, randomIA, ~] = unique(reshape(wspktrain,run.nPC,[])', 'rows');
            % indices to each of the unique states in my sample (all runs):
            voteState{1,stc} = zeros(m, size(clusterStatesAll{stc},1));
            for ru = 1:m-1
                ru
                    [~,locStates] = ismember(wspktrain(:,:,ru)', clusterStatesAll{stc}, 'rows');
                    % Gamhse...
                    voteState{1,stc}(ru,:) = accumarray([1:numel(voteState{1,stc}(ru,:)) locStates'].',[voteState{1,stc}(ru,:) ones(1,length(locStates))])' ;
            end
            [~,sidx] = sort(mean(voteState{1,stc}),'descend');
            figure;imagesc(voteState{1,stc}(:,sidx(2:100)))
            figure;plot(sum(voteState{1,stc}(:,sidx(2:100)))/size(clusterStatesAll{stc},1))
        end
        
        for k = 1:n_over_b
            ts(k) = ((k-1) * b) + 1 ;
            te(k) = (k) * b ;
            tmp_xi = sum(wspktrain(:,ts(k):te(k),:),2) ./ b ;
%             tmp_xi(tmp_xi==0) = 1e-12;
            bar_Xi(:,k,:) = tmp_xi;
            tmp_x = sum(sum(wspktrain(:,ts(k):te(k),:),2),3) ./ (b * m) ;
%             tmp_x(tmp_x==0) = 1e-11;
            dd_X(:,k) = tmp_x;
        end

        % gia SURF GRAPH!
        [~,IDX] = sort(mean(fr(labels==1,:),2),'descend')
        tmp = find(labels==1);
        cellpool = {};
        % prepei na elegxw oti isxyoun :
        % a) den exw 0 activity sto mean twn Q tou kane batch
        % b) den exw 0 activity sto parapanw, gia ka8e chain
        for kk = 2:2:length(tmp)
            tmp2 = tmp(IDX(1:kk));
            % a)
            tmp3 = ones(size(tmp2));
%             for ch=1:m
                for k = 1:n_over_b
                    tmp3 = tmp3 & any(wspktrain(tmp2,ts(k):te(k),ch), 2);
                end
%             end
            % b)
            if any(tmp3)
                tmp4 = tmp2(tmp3);
            else
                continue;
            end
            cellpool = [cellpool,tmp4];
        end
        mcs = cellfun(@(x)(mat2str(x)),cellpool,'uniformoutput',false);
        [un, idx2, idx] = unique(mcs);
        [~,sidx] = sort(cellfun(@length,cellpool(:,idx2)))
        cellpool = cellpool(:,idx2);
        cellpool = cellpool(:,sidx);

        for cp = 1:size(cellpool,2)
            fprintf('@ Q: %d, b: %d, CP: %d\n',Q, b, cp) ;
            % to xlen prepei na einai 700!! swsta?
            x_len = size(cellpool{1,cp},1) ;
%             % Calculate autocorrelation:
%             A = num2cell(wspktrain(cellpool{1,cp},:,1),1) ;
%             B = cell2mat(cellfun(@(x) bi2de(x'),A,'UniformOutput',false) ) ;
%             [acor,lag] = xcorr(B);
            % Multivariate extension as in: Brooks and Gelman, 1998:
            for k = 1:n_over_b
                ts(k)
%                 ts(k) = ((k-1) * b) + 1 ;
%                 te(k) = (k) * b ;
                cumSum_W = zeros(x_len,x_len,m) ;
                for ch = 1:m
                    cumSum_W(:,:,ch) = sum(reshape(cell2mat(cellfun(@(x) bsxfun(@times,x,x'),num2cell(bsxfun(@minus, wspktrain(cellpool{1,cp},ts(k):te(k),ch), bar_Xi(cellpool{1,cp},k,ch)),1),'uniformoutput', false)),x_len,x_len,[]),3) ;
                    bar_a{cp,k}(ch,:) = cellfun(@(x)bi2de(x'),num2cell(wspktrain(cellpool{1,cp},ts(k):te(k),ch),1),'uniformoutput',true );
                end
                % Instead of n we care about batch size (b):
                GR_W{cp,k} = sum(cumSum_W,3) ./ (m * (b-1));
                % reshape below might not be necessary..
%                 psicell = cellfun(@(x) bsxfun(@times,x,x'),num2cell(bsxfun(@minus, reshape(bar_Xi(cellpool{1,cp},k,:),x_len,[]), dd_X(cellpool{1,cp},k)),1),'uniformoutput', false);
%                 psi_d = reshape(bar_Xi(cellpool{1,cp},k,1),x_len,[]);
%                 psi_dd = dd_X(cellpool{1,cp},k);
%                 psi_minus = psi_d-psi_dd;
%                 psi_mult = psi_minus * psi_minus' ;
%                 psimat = cell2mat(psicell);
%                 psireshaped = reshape(psimat,x_len,x_len,[]);
% %                 % Memmory efficient algorithm:
% %                 psi_sub = num2cell(bsxfun(@minus, reshape(bar_Xi(cellpool{1,cp},k,:),x_len,[]), dd_X(cellpool{1,cp},k)),1);
% %                 psi_sub_mat=zeros( size(cellpool{1,cp},1));
% %                 for ii = 1:size(psi_sub,2)
% %                     psi_sub_mat = psi_sub_mat + psi_sub{1,ii} * psi_sub{1,ii}';
% %                 end
% %                 GR_B_N{cp,k} =  ./ (m -1);
% %                 
% %                 
                GR_B_N{cp,k} = sum(reshape(cell2mat(cellfun(@(x) bsxfun(@times,x,x'),num2cell(bsxfun(@minus, reshape(bar_Xi(cellpool{1,cp},k,:),x_len,[]), dd_X(cellpool{1,cp},k)),1),'uniformoutput', false)),x_len,x_len,[]),3) ./ (m -1);
                hat_V{cp,k} = (((b-1)/b)*GR_W{cp,k}) + ((1+(1/m))*GR_B_N{cp,k});
                
                if rcond(GR_W{cp,k}) < 1e-12
                    fprintf('Damn it. W matrix is near singular @t=%d. Skipping this t...\n',ts(k))
%                     tmp_mat = GR_W{cp,ig}\GR_B_N{cp,ig};
%                     tmp_mat(isnan(tmp_mat)) = 0;
%                     try hat_R(cp,ig) = ((N-1)/N) + ((M+1)/M) * max(eig(tmp_mat));
%                     catch ex
                        hat_R(cp,k) = NaN;
%                     end                    
                    % Since determinant is not an option, use something
                    % else for distance between matrices:
%                     dist_WV(cp,k) = sum(sum(abs(GR_W{cp,k} -  hat_V{cp,k})));
                    continue;
                else
                    hat_R(cp,k) = ((b-1)/b) + ((m+1)/m) * max(eig(GR_W{cp,k}\GR_B_N{cp,k}));
                end
            end
            ALL_hat_R{Ns,Qs,cp} = hat_R(cp,:) ;
            ALL_GR_B_N{Ns,Qs,cp} = GR_B_N(cp,:) ;
            ALL_GR_W{Ns,Qs,cp} = GR_W(cp,:) ;
            ALL_hat_V{Ns,Qs,cp} = hat_V(cp,:) ;
            ALL_bar_a{Ns,Qs,cp} = bar_a(cp,:) ;
        end
    end
end
%% Normalize range!
tstop = 20000;
nrange = cell(length(Nseq),length(Qseq));
for Qs = 1:length(Qseq)
    for Ns = 1:length(Nseq)
        Q = Qseq(Qs) ; % simple window (ms)
        Qr = floor(n / Q) ; % length of wspiketrain
%         n_over_b = 10 ; % No of batches
        b = floor(Qr/n_over_b); % Batches' length (b)

        gstart = (Q*b)/2 ;
        gend = tstop - (Q*b) ;
        nrange{Qs,Ns} = linspace(gstart,gend,n_over_b);
    end
end
% for Qs = 1:length(Qseq)
%     for Ns = 1:length(Nseq)
%         Q = Qseq(Qs) ; % simple window (ms)
% %         N = Nseq(Ns) ; % group length (n) after applying the window!
%         N = round((tstop/Q)/20);
%         Qr = floor(run.tstop / Q) ; % length of wspiketrain
%         ng = floor(Qr / N) ; % No of batches
%         gstart = (Q*N)/2 ;
%         gend = tstop - (Q*N) ;
%         nrange{Qs,Ns} = linspace(gstart,gend,ng);
%     end
% end


%Different Qs
figure; hold on;
cm = hsv(length(Qseq));
Ns = 1
lc = cell(1,length(Qseq));
for Qs = 1:length(Qseq)
    plot(nrange{Qs,Ns},ALL_hat_R{Ns,Qs,1},'color',cm(Qs,1:3));
    lc(1,Qs) = {num2str(Qseq(Qs))};
end
title(sprintf('N : %d batches',20))
legend(lc);


figure; hold on;
cm = hsv(length(Qseq));
Ns = 1
lc = cell(1,length(Qseq));
for Qs = 1:length(Qseq)
    myvar = [];
    for ii = 1:size(ALL_GR_B_N{Ns,Qs,cp},2)
        myvar(ii) = sum(sum(abs(ALL_GR_B_N{Ns,Qs,cp}{ii} - ALL_GR_W{Ns,Qs,cp}{ii})));
    end
    plot(nrange{Qs,Ns},myvar,'color',cm(Qs,1:3));
    lc(1,Qs) = {num2str(Qseq(Qs))};
end
title(sprintf('N : %d batches',20))
legend(lc);

% % make zero values black:
% var = (max(ALL_GR_B_N{Ns,Qs,1}{1}(:)) - min(ALL_GR_B_N{Ns,Qs,1}{1}(:)))/1000;
% var2 = min(ALL_GR_B_N{Ns,Qs,1}{1}(:)):var:max(ALL_GR_B_N{Ns,Qs,1}{1}(:));
% cm = hsv( length(unique(ALL_GR_B_N{Ns,Qs,1}{1}(:))) );
% figure;imagesc(ALL_GR_B_N{Ns,Qs,1}{1})
% colormap(cm)

figure; hold on;
cm = hsv(length(Qseq));
Ns = 1
lc = cell(1,length(Qseq));
for Qs = 1:length(Qseq)
    plot(nrange{Qs,Ns},cellfun(@(x) det(x),ALL_GR_B_N{Ns,Qs,1}),'color',cm(Qs,1:3));
    lc(1,Qs) = {num2str(Qseq(Qs))};
end
title(sprintf('N : %d batches',20))
legend(lc);

%Different Ns
figure; hold on;
cm = hsv(length(Nseq));
Qs = 1
lc = cell(1,length(Nseq));
for Ns = 1:length(Nseq)
    plot(nrange{Qs,Ns},ALL_hat_R{Ns,Qs,1},'color',cm(Ns,1:3));
    lc(1,Ns) = {num2str(Nseq(Ns))};
end
title(sprintf('Q : %dms',Qseq(Qs)))
legend(lc);
%% DEBUG wspktrain
wspktrain(run,Qseq,Nseq,Bin,stc,RUNS)
%% Compute ACT (autocorrelation time) as in :
% A Comparison of Methods for Computing Autocorrelation Time
% Also check Fishman, 1987.
% Effective Sample Size (ESS) and mixing measure are inverse proportional
% to the ACT.

M = run.nruns ; %No of chains (m)
Nseq = 10:50:700 ;
Qseq = 2:2:20;
Bin = 0.05:0.05:0.5 ; % Percentage of burn in samples.
ACT = zeros(length(Bin),length(Nseq),length(Qseq),size(cellpool,2),M) ;
Alpha = cell(length(Bin),length(Nseq),length(Qseq),size(cellpool,2),M) ;

for Bs = 1:length(Bin)
    for Qs = 1:length(Qseq)
        for Ns = 1:length(Nseq)
            Q = Qseq(Qs) ; % simple window (ms)
            N = Nseq(Ns) ; % group length (n) after applying the window!
            Bi = Bin(Bs) ;
            bnstart = round(run.tstop * Bi);
            brnrlen = run.tstop - bnstart ; % Burned run length
            Qr = floor(brnrlen / Q) ; % length of wspiketrain
            ng = floor(Qr / N) ; % No of groups
            bar_a = cell(size(cellpool,2),ng);
            s = (((1:Qr)-1) * Q) + 1 ;
            e = ((1:Qr) * Q) ;
            ts = [];
            te = [];
            
            tic;
            % Calculate spike trains ONLY for the above selected cells:
            
            wspktrain = zeros(run.nPC,Qr,M);
            for ru = 1:M
                spktrain = zeros(run.nPC,brnrlen);
                for c=1:run.nPC 
                    spks = round(RUNS{1,stc}{c,ru}.spikes') ;
                    spks = spks(spks > bnstart) - bnstart ;
                    spktrain(c,spks) = 1;
                end
                for k=1:Qr
                    wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
                end
            end
            runtime = toc
            
            for ig = 1:ng
                ts(ig) = ((ig-1) * N) + 1 ;
                te(ig) = (ig) * N ;
            end
            
            cells = find(labels==stc)';
            cellpool{1,1} = cells;
            for cp = 1:1%size(cellpool,2)
                % Multivariate extension as in: Brooks and Gelman, 1998:
                for ig = 1:ng
                    ts(ig) = ((ig-1) * N) + 1 ;
                    te(ig) = (ig) * N ;
                    for ch = 1:M
                        %                     alpha_i{cp,ig,ch} = wspktrain(cellpool{1,cp},ts(ig):te(ig),ch) ;
                        %                     bar_alpha_i{cp,ig,ch} = mean(wspktrain(cellpool{1,cp},ts(ig):te(ig),ch),2) ;
                        bar_a{ig,ch} = cellfun(@(x)bi2de(x'),num2cell(wspktrain(cellpool{1,cp},ts(ig):te(ig),ch),1),'uniformoutput',true );
                    end
                end
                for ch = 1:M
                    Alpha{Bs,Ns,Qs,cp,ch} = reshape(cell2mat(bar_a(:,ch))',1,[]) ;
                end
                for ch = 1:M
                    s_squr_m = var(cellfun(@mean,bar_a(:,ch) )) ; % variance of batch means
                    s_squr = var(cell2mat( bar_a(:,ch)' )) ;
                    ACT(Bs,Ns,Qs,cp,ch) = ng * (s_squr_m / s_squr) ;
                end
            end
            fprintf('@ Brn: %d, Q: %d, N: %d, CP: %d Runtime: %f\n',Bi * 100, Q, N, cp, toc) ;
        end
    end
end

%% Plot ACT
% tmp=[];
% for k=1:length(Bin)
% tmp(:,:,k) = mean(reshape(ACT(k,:,:,1,:),length(Nseq),length(Qseq),[]),3) ;
% imagesc( tmp(2:end,:,k) );
% set(gca,'YDir','normal');
% pause(0.5);
% cla;
% end

mintmp=[];
figure();hold on;
for k=1:length(Nseq)
    for l=1:length(Qseq)
%         plot( (shiftdim(tmp(k,l,: ))) );
        [~,mintmp(k,l)] = min(shiftdim(tmp(k,l,: ))) ;
    end
end
imagesc(mintmp);
set(gca,'YDir','normal');
%% Measure CV
for c=1:run.nPC
var = diff(RUNS_str{1,stc}{c,1}.spikes);
    RUN_CV(c) = std(var)/mean(var);
end
figure;bar(0:0.1:2,histc(RUN_CV,0:0.1:2))

%%
% Edw to afisa:@ Q: 10, N: 560, CP: 3
M = run.nruns ; %No of chains (m)
Nseq = [10,20,50,100,200,500,1000];%10:50:700 ;
Qseq = [2,5,10,20];%2:5:20;

Qseq = 5:5:40;
% Nseq = round((tstop/Qseq)/4);
Nseq = 1;
ALL_hat_R = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_GR_B_N = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_GR_W = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_hat_V = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;
ALL_bar_a = cell(length(Nseq),length(Qseq),size(cellpool,2)) ;


for Qs = 1:length(Qseq) 
    for Ns = 1:length(Nseq) 
        Q = Qseq(Qs) ; % simple window (ms)
%         N = Nseq(Ns) ; % group length (n) after applying the window!
        % N is the size of batch in Qs and NOT the nubmer of batches!
        N = round((tstop/Q)/20);
        Qr = floor(run.tstop / Q) ; % length of wspiketrain
        ng = floor(Qr / N) ; % No of batches
        hat_R = zeros(size(cellpool,2),ng) ;
        GR_W = cell(size(cellpool,2),ng) ;
        GR_B_N = cell(size(cellpool,2),ng) ;
        hat_V = cell(size(cellpool,2),ng) ;
        bar_Xi = zeros(run.nPC,ng,M) ;
        dd_X = zeros(run.nPC,ng) ;
        bar_a = cell(size(cellpool,2),ng);
        s = (((1:Qr)-1) * Q) + 1 ;
        e = ((1:Qr) * Q) ;
        ts = [];
        te = [];
        x_len = [];
        
        % Calculate spike trains ONLY for the above selected cells:
        wspktrain = zeros(run.nPC,Qr,M);
        for ru = 1:M
            spktrain = zeros(run.nPC,run.tstop);
            for c=1:run.nPC % must be row vector!!
                if ~isempty(RUNS{1,stc}{c,ru})
                    spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
                end
            end
            for k=1:Qr
                wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
            end
        end        
%         A = num2cell(wspktrain(,1:1000,1),1) ;
%         B = cell2mat(cellfun(@(x) bi2de(x'),A,'UniformOutput',false) ) ;
%         
%         [acor,lag] = xcorr(B);
%         plot(lag,acor) ;
        
        for ig = 1:ng
            ts(ig) = ((ig-1) * N) + 1 ;
            te(ig) = (ig) * N ;
            tmp_xi = sum(wspktrain(:,ts(ig):te(ig),:),2) ./ N ;
            tmp_xi(tmp_xi==0) = 1e-12;
            bar_Xi(:,ig,:) = tmp_xi;
            tmp_x = sum(sum(wspktrain(:,ts(ig):te(ig),:),2),3) ./ (N * M) ;
            tmp_x(tmp_x==0) = 1e-11;
            dd_X(:,ig) = tmp_x;
        end
%         cells = find(labels==stc)';
%             cellpool{1,1} = cells;
        for cp = 1%:size(cellpool,2)
            fprintf('@ Q: %d, N: %d, CP: %d\n',Q, N, cp) ;
            x_len = size(cellpool{1,cp},2) ;
%             % Calculate autocorrelation:
%             A = num2cell(wspktrain(cellpool{1,cp},:,1),1) ;
%             B = cell2mat(cellfun(@(x) bi2de(x'),A,'UniformOutput',false) ) ;
%             [acor,lag] = xcorr(B);
            % Multivariate extension as in: Brooks and Gelman, 1998:
            for ig = 1:ng
                ts(ig)
                ts(ig) = ((ig-1) * N) + 1 ;
                te(ig) = (ig) * N ;
                cumSum_W = zeros(x_len,x_len,M) ;
                for ch = 1:M
                    cumSum_W(:,:,ch) = sum(reshape(cell2mat(cellfun(@(x) bsxfun(@times,x,x'),num2cell(bsxfun(@minus, wspktrain(cellpool{1,cp},ts(ig):te(ig),ch), bar_Xi(cellpool{1,cp},ig,ch)),1),'uniformoutput', false)),x_len,x_len,[]),3) ;
                    bar_a{cp,ig}(ch,:) = cellfun(@(x)bi2de(x'),num2cell(wspktrain(cellpool{1,cp},ts(ig):te(ig),ch),1),'uniformoutput',true );
                end
                GR_W{cp,ig} = sum(cumSum_W,3) ./ (M * (N-1));
                % reshape below might not be necessary..
                GR_B_N{cp,ig} = sum(reshape(cell2mat(cellfun(@(x) bsxfun(@times,x,x'),num2cell(bsxfun(@minus, reshape(bar_Xi(cellpool{1,cp},ig,:),x_len,[]), dd_X(cellpool{1,cp},ig)),1),'uniformoutput', false)),x_len,x_len,[]),3) ./ (M -1);
                hat_V{cp,ig} = (((N-1)/N)*GR_W{cp,ig}) + ((1+(1/M))*GR_B_N{cp,ig});
                
                if rcond(GR_W{cp,ig}) < 1e-12
                    fprintf('Damn it. W matrix is near singular @t=%d. Skipping this t...\n',ts(ig))
%                     tmp_mat = GR_W{cp,ig}\GR_B_N{cp,ig};
%                     tmp_mat(isnan(tmp_mat)) = 0;
%                     try hat_R(cp,ig) = ((N-1)/N) + ((M+1)/M) * max(eig(tmp_mat));
%                     catch ex
                        hat_R(cp,ig) = NaN;
%                     end                    
                    continue;
                else
                    hat_R(cp,ig) = ((N-1)/N) + ((M+1)/M) * max(eig(GR_W{cp,ig}\GR_B_N{cp,ig}));
                end
            end
            ALL_hat_R{Ns,Qs,cp} = hat_R(cp,:) ;
            ALL_GR_B_N{Ns,Qs,cp} = GR_B_N(cp,:) ;
            ALL_GR_W{Ns,Qs,cp} = GR_W(cp,:) ;
            ALL_hat_V{Ns,Qs,cp} = hat_V(cp,:) ;
            ALL_bar_a{Ns,Qs,cp} = bar_a(cp,:) ;
        end
    end
end
plot(hat_R(cp,:),'r')

%%
pathprefix = 'X:\Documents\Glia\';
save(sprintf('%s/RND_MSRF_WorkingLabels_700.mat',pathprefix),'ALL_hat_R','ALL_GR_B_N','ALL_GR_W','ALL_hat_V','ALL_bar_a','-v7.3');

%%
figure;hold on;
for k = 1:length(Qseq)
    plot(ALL_hat_R{1 ,k,1});
end
ALL_hat_R{1:length(Nseq) ,1:length(Qseq),cp}


%% Autocorrelation:
r = xcorr(wspktrain(:,:,1))
%%
close all;
fh = figure('NumberTitle','off','Name',sprintf('Stimulated Cluster: %d, Q = %d, N = %d, Batches = %d',stc,Q,N, ng));
fill([0,0,ng*Q*N,ng*Q*N],[1,1.1,1.1,1],[0.95,0.95,0.95],'edgecolor',[0.95,0.95,0.95]) ;hold on;

sih = plot(((1:ng) * Q*N) + 1-(Q*N/2) ,hat_R(stc,:),'linewidth',3,'color','k');hold on;
ih = plot(((1:ng) * Q*N) + 1-(Q*N/2) ,hat_R([1:stc-1,stc+1:end],:),'linewidth',1);hold on;
axis([0,run.tstop,min(hat_R(:)),max(hat_R(:))]);

cn = cell(1,size(cellpool,2)-1)
for k =3:size(cellpool,2)+1
    cn{k} = sprintf('Cluster %d',k-2);
end
cn{1} = 'Convergence limits';
cn{2} = sprintf('* Cluster %d',stc);
cn(stc+2) = [];

legend(cn) ;
% saveas(fh,sprintf('%s/%s/STR_stc_%d_Q_%d_N_%d',pathprefix,run.path,stc,Q,N),'fig')
% legend(sih, sprintf('* Cluster %d',stc)) ; hold on;

%%
close all;
andplot = 1 ;
burnin = 5;
for Qs = 1:length(Qseq)
    for Ns = 1:length(Nseq)
        for cp = 1:size(cellpool,2)
            Q = Qseq(Qs) ; % simple window (ms)
            N = Nseq(Ns) ; % group length (n) after applying the window!
            Qr = floor(run.tstop / Q) ; % length of wspiketrain
            ng = floor(Qr / N) ; % No of groups
            if andplot
                fh = figure('NumberTitle','off','Name',sprintf('Stimulated Cluster: %d, Q = %d, N = %d, Batches = %d',stc,Q,N, ng));hold on;
                fill([0,0,ng*Q*N,ng*Q*N],[1,1.1,1.1,1],[0.95,0.95,0.95],'edgecolor',[0.95,0.95,0.95]) ;hold on;
            end
            
            X = [((1:ng) * Q*N) + 1-(Q*N/2)]' ;
            
            Rhat = ALL_hat_R{Ns,Qs,cp}(stc,:)' ;
            
%             if find(Rhat<1.3,1) == min( find(Rhat(burnin:end)<1.3,1), find(X(burnin:end) > run.stimend,1) )
            if (find(Rhat(burnin:end)<1.3,1) * N )  > run.stimend 
                start = (Rhat<1.3) ;
            else
                start = (X > run.stimend) ;
            end
            start(1:burnin-1) = 0;
            IDX = ( (~isnan(Rhat)) & start ) ;
            
            if sum(IDX) < 4
                its(Ns,Qs,cp) = NaN ;
                continue;
            end
            
            ft = fit(X(IDX) ,Rhat(IDX),'exp2') ;
            fs = min(X(IDX)) ;
            fe = max(X(IDX)) ;
            tmpits = fs + round(find( (ft(fs:0.01:fe) - 1.1 ) < eps,1) / 100) ;
            if ~isempty(tmpits)
                its(Ns,Qs,cp) = tmpits ;
            else
                its(Ns,Qs,cp) = NaN ;
            end
            if andplot
                scatter(X(IDX) ,Rhat(IDX),'r') ;
                scatter(X(~IDX) ,Rhat(~IDX),'b') ;
                scatter(its(Ns,Qs,cp),1.1,'*k') ;
                fplot(ft,[fs,fe],'k') ;
                axis([0,run.tstop,min(Rhat),max(Rhat)]) ;
                pause() ;
                close;
            end
        end
    end
end

for cp = 1:size(cellpool,2)
    imagesc(Qseq,Nseq,its(:,:,cp)) ;
    set(gca,'Ydir','normal') ;
%     colormap(gray);
    pause();
    cla;
end
%% Plot Normalized W, V, B/n, mean, variance
close all;

for k =1:size(cellpool,2)
    figure();
    tmp = cellfun(@(x)det(x),GR_W(k,:),'uniformoutput',true) ; 
    plot(((1:ng) * Q*N) + 1-(Q*N/2),tmp,'color','k');hold on;
    tmp = cellfun(@(x)det(x),hat_V(k,:),'uniformoutput',true) ;
    plot(((1:ng) * Q*N) + 1-(Q*N/2),tmp,'color','r');hold on;
    legend({'W','V'},'Interpreter','latex');
    figure();
    tmp = cellfun(@(x)det(x),GR_B_N(k,:),'uniformoutput',true) ;
    plot(((1:ng) * Q*N) + 1-(Q*N/2),tmp,'color','g');hold on;
    legend({'B/n'},'Interpreter','latex');
    figure();
    tmp = cellfun(@(x)mean(x(:)),bar_a(k,:),'uniformoutput',true) ;
    plot(((1:ng) * Q*N) + 1-(Q*N/2),tmp,'color','b');hold on;
    tmp = cellfun(@(x)var(x(:))/N,bar_a(k,:),'uniformoutput',true) ;
    plot(((1:ng) * Q*N) + 1-(Q*N/2),tmp,'color','c');hold on;
    legend({'\={a}', 'Var[\={a}]'},'Interpreter','latex');
    pause();
    close all;
end
%% Marginal probability of each cluster:
for k =1:size(cellpool,2)
    for ig = 1:ng
        ts = ((ig-1) * N) + 1 ;
        te = (ig) * N ;
%         plot(reshape(sum(wspktrain{1,k}(:,ts:te,:),2) / N,size(wspktrain{1,k},1),[]) ;
        mp(:,ig,:) = sum(wspktrain{1,k}(:,ts:te,:),2) / N ;
    end
    cm = jet(size(wspktrain{1,k},1)) ;
    for c = 1:size(wspktrain{1,k},1)
        for ru = 1:run.nruns
            plot(((1:ng) * Q*N) + 1-(Q*N/2),mp(13,:,ru)),'color',cm(c,1:3) ; hold on;
        end
    end
end




%%
close all;
Sid = 1;
MeanFFperCluster =[];


for ru = 1:run.nruns
    FFperCluster = NaN(run.NC_str,run.NC_str,run.nPC);
    for clu = 1:run.NC_str
        for c=1:run.nPC
            FFperCluster(RUNS_str{1,clu}{c,ru}.clusterID,clu,c) = RUNS_str{1,clu}{c,ru}.spike_count(run.stimend,run.tstop) / ((run.tstop-run.stimend)/1000) ;
        end
    end
    
    MeanFFperCluster(:,:,ru) = nanmean(FFperCluster,3);
end


for clu = 1:run.NC_str
    figure();
    bar3(reshape(MeanFFperCluster(:,clu,:),run.NC_str,run.nruns));
    axis([0,run.nruns+1,0,run.NC_str+1,0,150])
    view([270,0]);
end

%% Probability of recruitment:
close all;
PAP_recruit_mean_str = [];
PAP_recruit_mean_rnd = [];

if (exist('RUNS_str','var'))
    PAPerClust_str=zeros(run.NC_str(Sid),run.NC_str(Sid),run.nruns);
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_str(Sid)
            for clu=1:run.NC_str(Sid)
                for c=1:run.nPC
                    if(RUNS_str{1,stc}{c,ru}.clusterID == clu)
                        if(RUNS_str{1,stc}{c,ru}.persistentActivity)
                            PAPerClust_str(clu,stc,ru) = PAPerClust_str(clu,stc,ru)+1;
                        end
                    end
                end
            end
        end
    end
    PAP_recruit_str  = PAPerClust_str ./ repmat(run.cellsPerCluster_str',[1,run.NC_str,run.nruns]);
    PAP_recruit_mean_str = mean(PAP_recruit_str,3);
end

if (exist('RUNS_rnd','var'))
    PAPerClust_rnd=zeros(run.NC_rnd(Sid),run.NC_rnd(Sid),run.nruns);
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_rnd(Sid)
            for clu=1:run.NC_rnd(Sid)
                for c=1:run.nPC
                    if(RUNS_rnd{1,stc}{c,ru}.clusterID == clu)
                        if(RUNS_rnd{1,stc}{c,ru}.persistentActivity)
                            PAPerClust_rnd(clu,stc,ru) = PAPerClust_rnd(clu,stc,ru) +1;
                        end
                    end
                end
            end
        end
    end
    PAP_recruit_rnd  = PAPerClust_rnd ./ repmat(run.cellsPerCluster_rnd',[1,run.NC_rnd,run.nruns]);
    PAP_recruit_mean_rnd = mean(PAP_recruit_rnd,3);
end

maxvalue = 0;
maxvalue(logical(~isempty(PAP_recruit_mean_str))) = max([PAP_recruit_mean_str(:);maxvalue]);
maxvalue(logical(~isempty(PAP_recruit_mean_rnd))) = max([PAP_recruit_mean_rnd(:);maxvalue]);
ccm = jet( ceil(maxvalue*100)  ) ;

if (~isempty(PAP_recruit_mean_str))
    ccm(1,:) = [0,0,0];
    figure('Name','Recruitement Structured' );
    imagesc(PAP_recruit_mean_str);
    colormap(ccm(round(min(PAP_recruit_mean_str(:))*100)+1:round(max(PAP_recruit_mean_str(:))*100),:));
    colorbar;
    set(gcf,'Color',[1,1,1]);
end

if (~isempty(PAP_recruit_mean_rnd))
    ccm(1,:) = [0,0,0];
    figure('Name','Recruitement Random' );
    imagesc(PAP_recruit_mean_rnd);
    colormap(ccm(round(min(PAP_recruit_mean_rnd(:))*100)+1:round(max(PAP_recruit_mean_rnd(:))*100),:));
    colorbar;
    set(gcf,'Color',[1,1,1]);
end

%% PER NETWORK whole net Hamming distance: (Structured)
% close all;
% if enabled the Hamming distance contains all the cells (not the
% thresholded activation of each cluster) the threshold is defined in the
% next variable:
allCells = 0 ;
activeClusterThreshold = 0.5 ;
activationThreshold = 1;

PAPvect_str = zeros(run.nPC,run.NC_str,run.nruns);

if (allCells)
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_str(Sid)
            for c=1:run.nPC
                PAPvect_str(c,stc,ru) = RUNS_str{1,stc}{c,ru}.persistentActivity ;
            end
        end
    end
    
else
    PAPerClust_str=zeros(run.NC_str(Sid),run.NC_str(Sid),run.nruns);
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_str(Sid) % stimulated cluster
            for clu =1:run.NC_str(Sid) % cluster in net
                for c=1:run.nPC
                    if( (RUNS_str{1,stc}{c,ru}.clusterID == clu) && (RUNS_str{1,stc}{c,ru}.persistentActivity) )
                        PAPerClust_str(clu,stc,ru) = PAPerClust_str(clu,stc,ru)+1;
                    end
                end
            end
        end
    end
    PAPvect_str = (PAPerClust_str ./ repmat(run.cellsPerCluster_str',[1,run.NC_str,run.nruns])) >= activeClusterThreshold ;
end

% Remove the network states that had no PA in any run:
eyeidx  = repmat(eye(run.NC_str),[1,1,run.nruns]);
PAPvect_str( logical(eyeidx(:)) ) = 0;
[~,c]=find(mean(PAPvect_str,3));
[u,ia,~] = unique(c);
ia(end+1) = length(c)+1; ia = ia(2:end)-ia(1:end-1);
ia = ia >= activationThreshold;
PAPvect_str = PAPvect_str(:,u(ia),:);

HD_intra_str = zeros( run.NC_str , ((run.nruns^2)-run.nruns)/2 ) * NaN;
for stc=1:size(PAPvect_str,2)
    binchar='';
    binchar = '';
    for ru = 1:run.nruns
        binarr = PAPvect_str(:,stc,ru)';
        binchar(ru,:) = char(cell2mat(arrayfun(@(x) dec2bin(x,1),binarr,'UniformOutput', false)));
    end
    HD_intra_str(stc,:) = pdist(double(binchar),'minkowski',1);
end

% HD_intra_str(HD_intra_str~=0) = NaN;
MeanHD_intra_str = HD_intra_str(~isnan(HD_intra_str));
% MeanHD_intra_str = nanmean(HD_intra_str,2);


combs_clu_str = combntns(1:size(PAPvect_str,2),2);
combs_run = combntns(1:run.nruns,2);


HD_inter_str = zeros(size(combs_clu_str,1),size(combs_run,1)) * NaN;
for stc=1:size(combs_clu_str,1)
    binarr_a='';
    binarr_b='';
    binchar_a = '';
    binchar_b = '';
    for ru = 1:length(combs_run)
        binarr_a = PAPvect_str(:,combs_clu_str(stc,1),combs_run(ru,1))';
        binarr_b = PAPvect_str(:,combs_clu_str(stc,2),combs_run(ru,2))';
        
        binchar_a(ru,:) = char(cell2mat(arrayfun(@(x) dec2bin(x,1),binarr_a,'UniformOutput', false)));
        binchar_b(ru,:) = char(cell2mat(arrayfun(@(x) dec2bin(x,1),binarr_b,'UniformOutput', false)));
        
        HD_inter_str(stc,ru) = pdist(double([binchar_a(ru,:);binchar_b(ru,:)]),'minkowski',1);
    end    
end

% HD_inter_str(HD_inter_str==0) = NaN;
MeanHD_inter_str = HD_inter_str(~isnan(HD_inter_str));
% MeanHD_inter_str = nanmean(HD_inter_str,2);

AllValues_str = [MeanHD_intra_str(:);MeanHD_inter_str(:)];
ValuesGroups_str = [ones(length(MeanHD_intra_str),1);ones(length(MeanHD_inter_str),1)*2];

[p_val_str,ANOVA_table_str,mult_comp_stats_str] = anova1(AllValues_str,ValuesGroups_str)


%% PER NETWORK whole net Hamming distance: (Random)
close all;
% if enabled the Hamming distance contains all the cells (not the
% thresholded activation of each cluster) the threshold is defined in the
% next variable:
allCells = 1 ;
activeClusterThreshold = 0.5 ;
activationThreshold = 1;

PAPvect_rnd = zeros(run.nPC,run.NC_rnd,run.nruns);

if (allCells)
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_rnd(Sid)
            for c=1:run.nPC
                PAPvect_rnd(c,stc,ru) = RUNS_rnd{1,stc}{c,ru}.persistentActivity ;
            end
        end
    end
    
else
    PAPerClust_rnd = zeros(run.NC_rnd(Sid),run.NC_rnd(Sid),run.nruns);
    for ru = 1:run.nruns % what run of this cluster
        for stc =1:run.NC_rnd(Sid) % stimulated cluster
            for clu =1:run.NC_rnd(Sid) % cluster in net
                for c=1:run.nPC
                    if( (RUNS_rnd{1,stc}{c,ru}.clusterID == clu) && (RUNS_rnd{1,stc}{c,ru}.persistentActivity) )
                        PAPerClust_rnd(clu,stc,ru) = PAPerClust_rnd(clu,stc,ru)+1;
                    end
                end
            end
        end
    end
    PAPvect_rnd = (PAPerClust_rnd ./ repmat(run.cellsPerCluster_rnd',[1,run.NC_rnd,run.nruns])) >= activeClusterThreshold ;
end

% Remove the network states that had no PA in any run:
eyeidx  = repmat(eye(run.NC_rnd),[1,1,run.nruns]);
PAPvect_rnd( logical(eyeidx(:)) ) = 0;
[~,c]=find(mean(PAPvect_rnd,3));
[u,ia,~] = unique(c);
ia(end+1) = length(c)+1; ia = ia(2:end)-ia(1:end-1);
ia = ia >= activationThreshold;
PAPvect_rnd = PAPvect_rnd(:,u(ia),:);

HD_intra_rnd = zeros( run.NC_rnd , ((run.nruns^2)-run.nruns)/2 ) * NaN;
for stc=1:size(PAPvect_rnd,2)
    binchar='';
    binchar = '';
    for ru = 1:run.nruns
        binarr = PAPvect_rnd(:,stc,ru)';
        binchar(ru,:) = char(cell2mat(arrayfun(@(x) dec2bin(x,1),binarr,'UniformOutput', false)));
    end
    HD_intra_rnd(stc,:) = pdist(double(binchar),'minkowski',1);
end

% HD_intra_str(HD_intra_str~=0) = NaN;
MeanHD_intra_rnd = HD_intra_rnd(~isnan(HD_intra_rnd));
% MeanHD_intra_str = nanmean(HD_intra_str,2);


combs_clu_rnd = combntns(1:size(PAPvect_rnd,2),2);
combs_run = combntns(1:run.nruns,2);


HD_inter_rnd = zeros(size(combs_clu_rnd,1),size(combs_run,1)) * NaN;
for stc=1:size(combs_clu_rnd,1)
    binarr_a='';
    binarr_b='';
    binchar_a = [];
    binchar_b = [];
    for ru = 1:length(combs_run)
        binarr_a = PAPvect_rnd(:,combs_clu_rnd(stc,1),combs_run(ru,1))';
        binarr_b = PAPvect_rnd(:,combs_clu_rnd(stc,2),combs_run(ru,2))';
        
        binchar_a(ru,:) = cell2mat(arrayfun(@(x) dec2bin(x,1),binarr_a,'UniformOutput', false));
        binchar_b(ru,:) = cell2mat(arrayfun(@(x) dec2bin(x,1),binarr_b,'UniformOutput', false));
        
        HD_inter_rnd(stc,ru) = pdist(double([binchar_a(ru,:);binchar_b(ru,:)]),'minkowski',1);
    end    
end

% HD_inter_str(HD_inter_str==0) = NaN;
MeanHD_inter_rnd = HD_inter_rnd(~isnan(HD_inter_rnd));
% MeanHD_inter_str = nanmean(HD_inter_str,2);

AllValues_rnd = [MeanHD_intra_rnd(:);MeanHD_inter_rnd(:)];
ValuesGroups_rnd = [ones(length(MeanHD_intra_rnd),1);ones(length(MeanHD_inter_rnd),1)*2];

[p_val_rnd,ANOVA_table_rnd,mult_comp_stats_rnd] = anova1(AllValues_rnd,ValuesGroups_rnd)

%% Find Hamming distance between runs:

% PAPvect_str = cell(run.NC_str(Sid),run.NC_str(Sid),run.nruns);
% for ru = 1:run.nruns % what run of this cluster
%     for stc =1:run.NC_str(Sid)
%         for clu =1:run.NC_str(Sid)
%             for c=1:run.nPC
%                 if(RUNS_str{1,clu}{c,ru}.clusterID == stc)
%                     PAPvect_str{clu,stc,ru} = [PAPvect_str{clu,stc,ru}; RUNS_str{1,clu}{c,ru}.persistentActivity];
%                 end
%             end
%         end
%     end
% end
% 
% 
% HD_str = cell(run.NC_str(Sid),run.NC_str(Sid));
% for stc =1:run.NC_str(Sid)
%     for clu=1:run.NC_str(Sid)
%         binchar='';
%         binarr = [];
%         for ru = 1:run.nruns
%             binarr = PAPvect_str{clu,stc,ru}';
%             binchar(ru,:) = cell2mat(arrayfun(@(x) dec2bin(x,1),binarr,'UniformOutput', false)); % ginetai kai pio apla..
%         end
%         tmp = pdist(binchar,'minkowski',1);
%         tmp(tmp==0) = NaN;
%         HD_str{clu,stc} = tmp;
%     end
% end
% 
% % Mean HD for each cluster:
% MeanHD_str = cellfun(@(x) nanmean(x),HD_str);
% 
% 
% PAPvect_rnd = cell(run.NC_rnd(Sid),run.NC_rnd(Sid),run.nruns);
% for ru = 1:run.nruns % what run of this cluster
%     for stc =1:run.NC_rnd(Sid)
%         for clu =1:run.NC_rnd(Sid)
%             for c=1:run.nPC
%                 if(RUNS_rnd{1,clu}{c,ru}.clusterID == stc)
%                     PAPvect_rnd{clu,stc,ru} = [PAPvect_rnd{clu,stc,ru}; RUNS_rnd{1,clu}{c,ru}.persistentActivity];
%                 end
%             end
%         end
%     end
% end
% 
% 
% HD_rnd = cell(run.NC_rnd(Sid),run.NC_rnd(Sid));
% for stc =1:run.NC_rnd(Sid)
%     for clu=1:run.NC_rnd(Sid)
%         binchar='';
%         binarr = [];
%         for ru = 1:run.nruns
%             binarr = PAPvect_rnd{clu,stc,ru}';
%             binchar(ru,:) = cell2mat(arrayfun(@(x) dec2bin(x,1),binarr,'UniformOutput', false)); % ginetai kai pio apla..
%         end
%         tmp = pdist(binchar,'minkowski',1);
%         tmp(tmp==0) = NaN;
%         HD_rnd{clu,stc} = tmp;
%     end
% end
% 
% % Mean HD for each cluster:
% MeanHD_rnd = cellfun(@(x) nanmean(x),HD_rnd);
% 
% 
% % Cluster the mean HD per state:
% T_str = cluster(linkage(pdist(MeanHD_str)),'cutoff',0.5);
% T_rnd = cluster(linkage(pdist(MeanHD_rnd)),'cutoff',0.5);
% 
% length(unique(T_str)) / length(T_str) 
% length(unique(T_rnd)) / length(T_rnd) 





% % OLD Find Hamming distance between runs:
% binchar='';
% HD_str = zeros( run.NC_str , ((run.nruns^2)-run.nruns)/2 );
% for cl=1:run.NC_str
%     for r = 1:run.nruns
%         binarr = PAP_diff_recruit_str(:,cl,r)';
%         binchar(r,:) = cell2mat(arrayfun(@(x) dec2bin(x,8),binarr,'UniformOutput', false));
%     end
%     HD_str(cl,:) = pdist(binchar,'minkowski',1);
% end
% 
% binchar='';
% HD_rnd = zeros( run.NC_rnd , ((run.nruns^2)-run.nruns)/2 );
% for cl=1:run.NC_rnd
%     for r = 1:run.nruns
%         binarr = PAP_diff_recruit_rnd(:,cl,r)';
%         binchar(r,:) = cell2mat(arrayfun(@(x) dec2bin(x,8),binarr,'UniformOutput', false));
%     end
%     HD_rnd(cl,:) = pdist(binchar,'minkowski',1);
% end
% 
% figure('Name','Hamming distance between runs');
% plot(sum(HD_str,2),'g'); hold on;
% plot(sum(HD_rnd,2),'r'); 
% legend('Structured', 'Random')
% 
% % My metric: mean probability of recruitement / Hamming distance ratio:
% figure('Name','my metric between runs');
% plot(mean(PAP_recruit_mean_str) ./ sum(HD_str,2)','g'); hold on;
% plot(mean(PAP_recruit_mean_rnd) ./ sum(HD_rnd,2)','r'); 
% legend('Structured', 'Random')



% ccmd = jet( ceil(max([PAP_diff_recruit_str(:);PAP_diff_recruit_rnd(:)])*100)  ) ;
% figure('Name','STD Diff Recruitement Random' );
% ccmd(1,:) = [0,0,0];
% for i=1:run.NC_rnd
%     diffrec(:,i) = std(reshape(PAP_diff_recruit_rnd(:,i,:),run.NC_rnd,[]),0,2);
% end
% imagesc(diffrec);
% colormap(ccmd(1:round(max(PAP_diff_recruit_rnd(:))*100),:));
% colorbar;
% 
% 
% ccmd(1,:) = [0,0,0];
% figure('Name','STD Diff Recruitement Structured' );
% for i=1:run.NC_str
%     diffrec(:,i) = std(reshape(PAP_diff_recruit_str(:,i,:),run.NC_str,[]),0,2);
% end
% imagesc(diffrec);
% colormap(ccmd(1:round(max(PAP_diff_recruit_str(:))*100),:));
% colorbar;



%%  Panos newer data!
% STRUCTURED NET
SPK_str = cell(run.NC_str(Sid),run.NC_str(Sid),run.nruns);
[SPK_str(:)] = deal({{1}});
counting_str=ones(run.NC_str(Sid),run.NC_str(Sid),run.nruns) ;
for ru=1:run.nruns
    for cl=1:run.NC_str(Sid)
        for clu=1:run.NC_str(Sid)
            for c=1:run.nPC
                if(RUNS_str{1,clu}{c,ru}.clusterID == cl)
                    SPK_str{cl,clu,ru}(counting_str(cl,clu,ru)) = {RUNS_str{1,clu}{c,ru}.spikes};
                    counting_str(cl,clu,ru) = counting_str(cl,clu,ru) + 1;
                end
            end
        end
    end
end

% RANDOM NET
SPK_rnd = cell(run.NC_rnd(Sid),run.NC_rnd(Sid),run.nruns);
[SPK_rnd(:)] = deal({{1}});
counting_rnd=ones(run.NC_rnd(Sid),run.NC_rnd(Sid),run.nruns) ;
for ru=1:run.nruns
    for cl=1:run.NC_rnd(Sid)
        for clu=1:run.NC_rnd(Sid)
            for c=1:run.nPC
                if(RUNS_rnd{1,clu}{c,ru}.clusterID == cl)
                    SPK_rnd{cl,clu,ru}(counting_rnd(cl,clu,ru)) = {RUNS_rnd{1,clu}{c,ru}.spikes};
                    counting_rnd(cl,clu,ru) = counting_rnd(cl,clu,ru) + 1;
                end
            end
        end
    end
end
%%
figure();
for c=1:75
    [~,B]=RUNS_rnd{1,clu}{c,ru}.spike_count(1,run.tstop);
    plot([(B./10)],[0;1000./diff(B./10)],'Color',rand(1,3));hold on;
    tmp = RUNS_rnd{1,clu}{c,ru}.mv;
    tmp = tmp(1:10:end);
    plot(tmp,'k');hold on;
    pause;
    cla
end

%% generate state data (connectivity, positions, clusters etc)
nPC = 75;
%  ---- Rest of the cells ----
nPV = round(nPC*13/75);
nCB = round(nPC*6/75);
nCR = round(nPC*6/75);
if(nPV==0)% NEURON ISSUE:
    nPV = 1;
end
if(nCB == 0)
    nCB = 1;
end
if(nCR == 0)
    nCR=1;
end
nAll = sum([nPC,nPV,nCB,nCR]);


connProbRange = 0.7;
totalStateNo = 10;

state_str = zeros(nPC,nPC,totalStateNo);
state_rnd = zeros(nPC,nPC,totalStateNo);

for sts = 1:totalStateNo
    gparams(sts).nPC = nPC;
    gparams(sts).nPV = nPV;
    gparams(sts).nCB = nCB;
    gparams(sts).nCR = nCR;
    gparams(sts).nAll = nAll;

    gparams(sts).connBinsPC2PC = 0:30:500;
    gparams(sts).connProbsPC2PC = 0.25 .* exp(-0.006 * gparams(sts).connBinsPC2PC);
    gparams(sts).connBinsPV2PC = 0:20:500;
    gparams(sts).connProbsPV2PC = linspace(1,0,length(gparams(sts).connBinsPV2PC)) * 0.3; % OVERRIDE!! Original was 0.7
    gparams(sts).connBinsCB2PC = 0:20:500;
    gparams(sts).connProbsCB2PC = linspace(1,0,length(gparams(sts).connBinsCB2PC)) ;
    gparams(sts).recipBinsPC2PC = 0:30:500;
    gparams(sts).recipProbsPC2PC = 0.12 .* exp(-0.006 * gparams(sts).recipBinsPC2PC);
    % Update probabilities from Perin et al. (Figure S4)
    gparams(sts).NconnBins = [0,1,2,3,4];
    gparams(sts).NincomingProbs = [0.1, 0.2, 0.25, 0.4, 0.52] ;
    gparams(sts).NoutgoingProbs = [0.12, 0.25, 0.24, 0.2, 0.23] ;
    
    
    %  ---- Initialize Pyramidal cells ----
    gparams(sts).PCsomata = CreateRandomNetwork(nPC, 90, 3); % was 50, to much...
    gparams(sts).distPC2PC = generateDistanceMat(gparams(sts).PCsomata', 0);
    
    % PV somata: about 40 per mm2 for 50um depth of slice as in:
    %      Ohira, K., Takeuchi, R., Iwanaga, T., & Miyakawa, T. (2013).
    %      Chronic fluoxetine treatment reduces parvalbumin expression and
    %      perineuronal nets in gamma-aminobutyric acidergic interneurons of
    %      the frontal cortex in adult mice. Molecular Brain, 6(1), 43.
    %      doi:10.1186/1756-6606-6-43
    
    % Ara exw peripou 12.5 PV se enan kybo 250x250x250 um
    gparams(sts).PVsomata = CreateCubeNetworkPV(90, nPV); % 226.6 per mm squared (?!?) paper???
    gparams(sts).CBsomata = CreateCubeNetworkPV(90, nCB);
    gparams(sts).CRsomata = CreateCubeNetworkPV(0, nCR);
    
    % Distance of each interneuron from Pyramidal and connect based on
    % probability from Yuste 2011:
    
    % gia to connectivity twn PV 2 PC  o Packer et al., 2011:
    % elegxei gia ena region gyrw apo to ka8e PC. Epomenws ta histograms einai
    % swsta gia ta PV2PC (ka8ws briskontai sto idio layer). Oso gia ta
    % CB(SOM)2PC pou briskontai se diaforetika layers mporoume na eikasoume
    % oti:
    % * Apo to sxhma twn SOM (Packer et al., 2013), oso metakineisai sto
    % transverse layer oi tuft dendrites enws PC 8a exoun tis idies pi8anotites
    % gia overlap opos exoun kai ta PC pou briskontai sto idio layer. Epomenws
    % metraw san distance CB2PC mono to distance sto transverse plane kai
    % xrisimopoiw tis pi8anotites pou dinei o Packer et al., 2013, Fig4D
    gparams(sts).distPV2PC = distancePV2PC(gparams(sts).PVsomata,gparams(sts).PCsomata);
    gparams(sts).distCB2PC = distancePV2PC(gparams(sts).CBsomata,gparams(sts).PCsomata);
    
    gparams(sts).ConnMatPV2PC = connectPV2PC(gparams(sts).distPV2PC,gparams(sts).connBinsPV2PC,gparams(sts).connProbsPV2PC);
    gparams(sts).ConnMatCB2PC = connectCB2PC(gparams(sts).distCB2PC,gparams(sts).connBinsCB2PC,gparams(sts).connProbsCB2PC);
    gparams(sts).CRsomata = CreateCubeNetworkPV(0, nCR);
    gparams(sts).ConnMatPC2PV = connectPC2PV(gparams(sts).distPV2PC')  ;
    %             constraintPVreciprocity();
    
    
    % Takeshi Otsuka, 2009:
    % Test if the PC2PV reciprocals are 20% of PC2PV pairs:
    ctr = 10000;
    while ctr && (0.1 < abs(((sum(sum((gparams(sts).ConnMatPC2PV & gparams(sts).ConnMatPV2PC')))*100) / numel(gparams(sts).ConnMatPC2PV)) - 20))
        % if reciprocals are more/less, move connections to reach 20%:
        foundIt = 0;
        while ~foundIt
            ip = ceil(rand(1)*size(gparams(sts).ConnMatPC2PV,1));
            jp = ceil(rand(1)*size(gparams(sts).ConnMatPC2PV,2));
            if gparams(sts).ConnMatPC2PV(ip,jp) && ~gparams(sts).ConnMatPV2PC(jp,ip)
                foundIt = 1;
            end
        end
        foundIt = 0;
        while ~foundIt
            ii = ceil(rand(1)*size(gparams(sts).ConnMatPC2PV,1));
            ji = ceil(rand(1)*size(gparams(sts).ConnMatPC2PV,2));
            if ~gparams(sts).ConnMatPC2PV(ii,ji) && gparams(sts).ConnMatPV2PC(ji,ii)
                foundIt = 1;
            end
        end
        % move the projection of the PC to form reciprocal with the PV:
        gparams(sts).ConnMatPC2PV(ip,jp) = 0;
        gparams(sts).ConnMatPC2PV(ii,ji) = 1;
        ctr = ctr - 1 ;
    end
    % Takeshi Otsuka, 2009:
    % Test if the PC2PV connections are ~2% :
    gparams(sts).ConnMatPC2PV = gparams(sts).ConnMatPC2PV';
    tmpD = find(xor(gparams(sts).ConnMatPC2PV,gparams(sts).ConnMatPV2PC)) ;
    gparams(sts).ConnMatPC2PV(tmpD) = 0;
    gparams(sts).ConnMatPV2PC(tmpD) = 1;
    tmpR = rand(length(tmpD),1) < 0.02;
    gparams(sts).ConnMatPC2PV(tmpR) = 1;
    gparams(sts).ConnMatPV2PC(tmpR) = 0;
    gparams(sts).ConnMatPC2PV = gparams(sts).ConnMatPC2PV';
    
    % Set indices for ease of mind:
    pc = size(gparams(sts).PCsomata,1);
    pv = size(gparams(sts).PCsomata,1) + size(gparams(sts).PVsomata,1);
    cb = size(gparams(sts).PCsomata,1) + size(gparams(sts).PVsomata,1) + size(gparams(sts).CBsomata,1);
    cr = size(gparams(sts).PCsomata,1) + size(gparams(sts).PVsomata,1) + size(gparams(sts).CBsomata,1) + size(gparams(sts).CRsomata,1);
    
    % Populate the final all-to-all connectivity matrix:
    gparams(sts).connmatrix = zeros(nAll);
    %         % Pyramidals to all types of interneurons:
    %         AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
    % PCs to PVs
    gparams(sts).connmatrix(1:pc,pc+1:pv) = gparams(sts).ConnMatPC2PV;
    %     AllConnMat(1:pc,pv+1:end) = 1; % Connected to all interneurons: too much. maybe less connectivity?
    gparams(sts).connmatrix(1:pc,pv+1:end) = rand(nPC, nAll-pv) > 0.9; % Connected to all interneurons: too much. maybe less connectivity?
    % PVs connect to all other PVs + autapses (need to gap junctions?)
    gparams(sts).connmatrix(pc+1:pv,pc+1:pv) = 1;
    % PVs to PCs based on above connectivity (Yuste 2011):
    gparams(sts).connmatrix(pc+1:pv,1:pc) = gparams(sts).ConnMatPV2PC;
    % CB only connect to PC (Xenia) na to psa3w...
    gparams(sts).connmatrix(pv+1:cb,1:pc) = gparams(sts).ConnMatCB2PC;
    % CR only connect to PC and CB (Xenia) na to psa3w...
    gparams(sts).connmatrix(cb+1:cr,1:pc) = 1;
    gparams(sts).connmatrix(cb+1:cr,pv+1:cb) = 1;
    
    
    
    factorOverallProb = connProbRange;
    % Pyramidals connect to all
    gparams(sts).connProbsPC2PC = gparams(sts).connProbsPC2PC*factorOverallProb;
    gparams(sts).recipProbsPC2PC = gparams(sts).recipProbsPC2PC *factorOverallProb;
    gparams(sts).NincomingProbs = gparams(sts).NincomingProbs *factorOverallProb;
    gparams(sts).NoutgoingProbs = gparams(sts).NoutgoingProbs *factorOverallProb;
    
    MAXPC2PC = [];
    CCDIFF = [];
    CCMAX = 0;
    tmpMat_str = [];
    tmpMat_rnd = [];
    for manytimes=1:200
        manytimes
        tmpMat_str = create_graph_CN(gparams(sts).distPC2PC,gparams(sts).connBinsPC2PC, gparams(sts).connProbsPC2PC, ...%0.07  %0.2 dinei PA
            gparams(sts).recipBinsPC2PC, gparams(sts).recipProbsPC2PC,gparams(sts).NconnBins,gparams(sts).NincomingProbs,gparams(sts).NoutgoingProbs);
        prob_conn_rnd_ind = sum(sum(tmpMat_str)) / (numel(tmpMat_str)-nPC); % joint probability
        tmpMat_rnd = create_graph_WS(gparams(sts).distPC2PC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
        [mcccs,~,cccs] = clust_coeff( tmpMat_str );
        [mcccr,~,cccr] = clust_coeff( tmpMat_rnd );
        CCDIFF(manytimes) = mcccs - mcccr;
        if(CCMAX < CCDIFF(manytimes))
            disp(sprintf('Replacing %f with %f',CCMAX,CCDIFF(manytimes)))
            CCMAX = CCDIFF(manytimes);
            MAXPC2PC(:,:,1) = tmpMat_str;
            MAXPC2PC(:,:,2) = tmpMat_rnd;
        end
    end
    max(CCDIFF)
    states.state_str(:,:,sts) = gparams(sts).connmatrix;
    states.state_rnd(:,:,sts) = gparams(sts).connmatrix;
    states.state_str(1:pc,1:pc,sts) = MAXPC2PC(:,:,1);
    states.state_rnd(1:pc,1:pc,sts) = MAXPC2PC(:,:,2);
    
    
    % Find nearest neighbors of Pyramidals only
    [CNi_str,CNo_str] = m_commonNeighbors(states.state_str(1:pc,1:pc,sts));
    [CNi_rnd,CNo_rnd] = m_commonNeighbors(states.state_rnd(1:pc,1:pc,sts));
    gparams(sts).mergedCN_str = [CNi_str + CNo_str] .* states.state_str(1:pc,1:pc,sts);
    gparams(sts).mergedCN_rnd = [CNi_rnd + CNo_rnd] .* states.state_rnd(1:pc,1:pc,sts);
    gparams(sts).mergedCN_strT = gparams(sts).mergedCN_str';
    gparams(sts).mergedCN_str(logical(tril(ones(size(gparams(sts).mergedCN_str)),-1))) = gparams(sts).mergedCN_strT(logical(tril(ones(size(gparams(sts).mergedCN_str)),-1)));
    gparams(sts).mergedCN_rndT = gparams(sts).mergedCN_rnd';
    gparams(sts).mergedCN_rnd(logical(tril(ones(size(gparams(sts).mergedCN_rnd)),-1))) = gparams(sts).mergedCN_strT(logical(tril(ones(size(gparams(sts).mergedCN_rnd)),-1)));
    
    % Affinity Propagation:
    gparams(sts).cellsPerCluster_str = [];
    gparams(sts).cellsPerCluster_rnd = [];
    
    %perform affinity propagation algorithm:
    Sid=1; % depricated variable from another clustering algo..
    
    % try to force different No of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    srtClNo = 10;
    setFlg = 1;
    while setFlg
        for i=1:10
            [idx_str,~,~,~,~]=apclusterK(gparams(sts).mergedCN_str,srtClNo);
            gparams(sts).NC_str(Sid) = length(unique(idx_str));
            [~,~,gparams(sts).labels_str(:,Sid)] = unique(idx_str);
            if(min(histc(gparams(sts).labels_str(:,Sid),1:gparams(sts).NC_str(Sid))') > 6)
                setFlg = 0;
                break;
            end
            
        end
        srtClNo = srtClNo - 1;
    end
    
    
    % try to force different nNo of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    srtClNo = 10;
    setFlg = 1;
    while setFlg
        for i=1:10
            [idx_rnd,~,~,~,~]=apclusterK(gparams(sts).mergedCN_rnd,srtClNo);
            gparams(sts).NC_rnd(Sid) = length(unique(idx_rnd));
            [~,~,gparams(sts).labels_rnd(:,Sid)] = unique(idx_rnd);
            if(min(histc(gparams(sts).labels_rnd(:,Sid),1:gparams(sts).NC_rnd(Sid))') > 6)
                setFlg = 0;
                break;
            end
            
        end
        srtClNo = srtClNo - 1;
    end
    
    
%     gparams(sts).NC_rnd(Sid) = length(unique(gparams(sts).labels_rnd(:,sts)));
%     [~,~,gparams(sts).labels_rnd(:,Sid)] = unique(gparams(sts).labels_rnd(:,sts));
%     obj.NC_str(Sid) = length(unique(states.Labels_str(:,obj.state)));
%     [~,~,obj.labels_str(:,Sid)] = unique(states.Labels_str(:,obj.state));
    %how many cells in each structured cluster?
    gparams(sts).cellsPerCluster_str = histc(gparams(sts).labels_str(:,Sid),1:gparams(sts).NC_str(Sid))';
    %how many cells in each random cluster?
    gparams(sts).cellsPerCluster_rnd = histc(gparams(sts).labels_rnd(:,Sid),1:gparams(sts).NC_rnd(Sid))';
    gparams(sts).stimCellsPerCluster = min(min(gparams(sts).cellsPerCluster_str),min(gparams(sts).cellsPerCluster_rnd));
            
end

states.gparams = gparams;
% Save them in a structure nice and tidy..

save('states_06_G.mat','states')