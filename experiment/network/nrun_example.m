clear all;close all;clc;
cd('C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network')
% cd('E:\NEURON_RUNS')
load('states_07_G.mat')
pathprefix = 'N:/NEURON_PROJECTS/NEW_RUNS/';

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
clustbiasRange = 0:0.1:1;
for sn = 1:10
    clearvars -except states sn tmpSPK pathprefix clustbiasRange
    % obj = nrun(experimentid,npyrs,serial,state,exprun,tstop)
    ID = 12;
    SN = sn;
    ST = 1;
    run = nrun(ID,75,SN,ST,10,5000);
    run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
    
    run.init(states);
    
    % Insert cynaptic clustering bias: 0 = clustered, 1 = random (no
    % clustering)
    run.CLUSTBIAS = 0;%clustbiasRange(sn);
    
%     sum(sum(run.stateSTR(1:run.nPC,1:run.nPC)))            / (run.nPC*run.nPC)
    save(sprintf('%s/EXP_ID%d_SN%d_ST%d.mat',run.path,run.id,run.sn,run.state),'run');
    save(sprintf('%s%s/EXP_ID%d_SN%d_ST%d.mat',pathprefix,run.path,run.id,run.sn,run.state),'run');
    
    
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
    run.exportParams2Neuron();

    % Push data to cluster:
    run.push();
    
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
load(sprintf('%sexperiment_%d/EXP_ID%d_SN%d_ST%d.mat',pathprefix,12,12,1,1)) % 8,10


RUNS_str = {};
Sid=1;
fprintf('Loading runs...');
PCcells_str = {};
for stc=1:run.NC_str(Sid)
    stc
    for ru = 1:run.nruns
        for c=1:run.nPC
            if (run.ISBINARY)
                mycell = ncell(nrn_vread(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,c-1,ru-1),'n'),10);
            else
                mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,stc-1,c-1,ru-1)),10);
            end
            %             mycell = ncell(load(sprintf('%sexperiment_10_10Hz_Stim/%s/%d_%d_%d.txt',mypath,'STR',t,c-1,ru-1)),10);
            mycell.clusterID = run.labels_str(c,Sid);
            %             mycell.position = PCsomata(c,1:3);
            %             PCcells_str{c,ru}=mycell.hasPersistent(run.stimend-1,25,run.tstop-(run.stimend-1)); % paper?
            [S,~,~] = findUPstates(mycell.mv(run.stimend*run.dt:run.dt:end),10, 3, -66, 3000 );
            if ~isempty(S);
                mycell.persistentActivity = 1;
            else
                mycell.persistentActivity = 0;
            end
            PCcells_str{c,ru} = mycell ;
            %             leastSynapses(c) = PCcells_str{c,ru}.nspikes;
        end
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
    %     all(leastSynapses)
end
fprintf('DONE!\n');
%% Load RND
close all; clear all; clc;
load(sprintf('STR_%d_%d.mat',12,2))

RUNS_rnd = {};
Sid=1;
fprintf('Loading runs...');
PCcells_rnd = {};
for clu=1:run.NC_rnd(Sid)
    clu
    for ru = 1:run.nruns
        for c=1:run.nPC
            if (run.ISBINARY)
                mycell = ncell(nrn_vread(sprintf('%s%s/RND_%d/%d_%d_%d.bin',pathprefix,run.path,run.state,clu-1,c-1,ru-1),'n'),10);
            else
                mycell = ncell(load(sprintf('%s/RND_%d/%d_%d_%d.txt',run.path,run.state,clu-1,c-1,ru-1)),10);
            end
            %             mycell = ncell(load(sprintf('%sexperiment_10_10Hz_Stim/%s/%d_%d_%d.txt',mypath,'STR',t,c-1,ru-1)),10);
            mycell.clusterID = run.labels_rnd(c,Sid);
            %             mycell.position = PCsomata(c,1:3);
            %             PCcells_str{c,ru}=mycell.hasPersistent(run.stimend-1,25,run.tstop-(run.stimend-1)); % paper?
            [S,~,~] = findUPstates(mycell.mv(run.stimend*run.dt:run.dt:end),10, 3, -66, 3000 );
            if ~isempty(S);
                mycell.persistentActivity = 1;
            else
                mycell.persistentActivity = 0;
            end
            PCcells_rnd{c,ru} = mycell ;
            %             leastSynapses(c) = PCcells_str{c,ru}.nspikes;
        end
        %         for c=run.nPC+1:run.nPC+run.nPV
        %             mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,t-1,c-1,ru-1)),10);
        %             PVcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        %         end
        %         for c=run.nPV+1:run.nPV+run.nCB
        %             mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,t-1,c-1,ru-1)),10);
        %             CBcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        %         end
        %         for c=run.nCB+1:run.nCB+run.nCR
        %             mycell = ncell(load(sprintf('%s/STR_%d/%d_%d_%d.txt',run.path,run.state,t-1,c-1,ru-1)),10);
        %             CRcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        %         end
    end
    RUNS_rnd{1,clu} = PCcells_rnd(:,:);
    %     all(leastSynapses)
end
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
ru =1 ;% what run of this cluster
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
                            PAPerClust_rnd(clu,stc,ru) = PAPerClust_rnd(clu,stc) +1;
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
allCells = 1 ;
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
allCells = 0 ;
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

save('states_07_G.mat','states')