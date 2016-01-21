
aa=10
% for aa = 5
    close all; clearvars -except aa; 
    randn('seed',aa);
%     addpath(genpath('~/Documents/MATLAB/'));
%     addpath(genpath('~/Documents/GitHub/prefrontal-micro/'));
%     addpath(genpath('~/Documents/GitHub/neurocommitto/'));
    mypath = '/srv/userdata/HOMES/stefanos/Desktop/prefrontal-micro/experiment/network/';
    mypath2 = '~/Desktop/prefrontal-micro/';
%     mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
%     mypath2 = '~/Documents/GitHub/prefrontal-micro/';

    
    exprun=1;
    stimend = 1500;
    nPC = 75;%225; %216;
%     clstr = 1; % Cluster type
    
% % Create tissue sanbox of given dimensions to fit all the pertitnent
% % neuronal types in it, in correct densities.
% tSandbox.side = 250;
% PVsomata = CreateCubeNetworkPV(250, 12.5);
%     scatter3(PVsomata(:,1),PVsomata(:,2),PVsomata(:,3));


    connBinsPC2PC = 0:30:500;
    connProbsPC2PC = 0.25 .* exp(-0.006 * connBinsPC2PC);
    connBinsPV2PC = 0:20:500;
%     connProbsPV2PC = linspace(1,0,length(connBinsPV2PC)) * 0.7; % wste na
%     mas volevei, POUSTIA....
    connProbsPV2PC = linspace(1,0,length(connBinsPV2PC)) * 0.7;
    connBinsCB2PC = 0:20:500;
%     connProbsCB2PC = normpdf(connBinsCB2PC,150,90) *10; %afto den einai
% connection probability, alla histogram of distances of connected CB2PC
% (layer2/3) pairs! (mporw na to xrisimopoiisw gia validation?? )
    connProbsCB2PC = linspace(1,0,length(connBinsCB2PC)) ;
    
    recipBinsPC2PC = 0:30:500;
    recipProbsPC2PC = 0.12 .* exp(-0.006 * recipBinsPC2PC);
    % Update probabilities from Perin et al. (Figure S4)
    NconnBins = [0,1,2,3,4];
    NincomingProbs = [0.1, 0.2, 0.25, 0.4, 0.52] ;
    NoutgoingProbs = [0.12, 0.25, 0.24, 0.2, 0.23] ;
    
    %  ---- Initialize Pyramidal cells ----
    PCsomata = CreateRandomNetwork(nPC, 90, 3); % was 50, to much...
    distPC2PC = generateDistanceMat(PCsomata', 0);
    scatter3(PCsomata(:,1),PCsomata(:,2),PCsomata(:,3));
    
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
    nAllCells = sum([nPC,nPV,nCB,nCR]);
    AllCells = [nPC , nPV , nCB , nCR, nAllCells];
    
      
    
    % PV somata: about 40 per mm2 for 50um depth of slice as in:
%      Ohira, K., Takeuchi, R., Iwanaga, T., & Miyakawa, T. (2013). 
%      Chronic fluoxetine treatment reduces parvalbumin expression and 
%      perineuronal nets in gamma-aminobutyric acidergic interneurons of 
%      the frontal cortex in adult mice. Molecular Brain, 6(1), 43. 
%      doi:10.1186/1756-6606-6-43
     
% Ara exw peripou 12.5 PV se enan kybo 250x250x250 um
    PVsomata = CreateCubeNetworkPV(90, nPV); % 226.6 per mm squared (?!?) paper???
    CBsomata = CreateCubeNetworkPV(90, nCB);
    CRsomata = CreateCubeNetworkPV(0, nCR);
    
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
    distPV2PC = distancePV2PC(PVsomata,PCsomata);
    distCB2PC = distancePV2PC(CBsomata,PCsomata);
    
    ConnMatPV2PC = connectPV2PC(distPV2PC,connBinsPV2PC,connProbsPV2PC);
    ConnMatCB2PC = connectCB2PC(distCB2PC,connBinsCB2PC,connProbsCB2PC);
    
       %     Packer, A. M., & Yuste, R. (2011). Dense, unspecific connectivity of 
% neocortical parvalbumin-positive interneurons: a canonical microcircuit 
% for inhibition? The Journal of Neuroscience : The Official Journal of the
% Society for Neuroscience, 31(37), 1326071. 
% doi:10.1523/JNEUROSCI.3131-11.2011

% Packer, A., McConnell, D., Fino, E., & Yuste, R. (2013). Axo-dendritic 
% overlap and laminar projection can explain interneuron connectivity to 
% pyramidal cells. Cerebral Cortex, (December), 27902802. 
% doi:10.1093/cercor/bhs210

    % Connection of FS interneurons to PCs as in :
    % Cortical inhibitory cell types differentially form intralaminar
    % and interlaminar subnetworks with excitatory neurons, Otsuka Takeshi
    % Kawaguchi Yasuo, 2009
    ConnMatPC2PV = connectPC2PV(distPV2PC')  ;
    % Test if the PC2PV reciprocals are 20% of PC2PV pairs:
    ctr = 10000;
    while ctr && (1.5 < abs(((sum(sum((ConnMatPC2PV & ConnMatPV2PC')))*100) / numel(ConnMatPC2PV)) - 20))
        % if reciprocals are more/less, move connections to reach 20%:
        foundIt = 0;
        while ~foundIt
            ip = ceil(rand(1)*size(ConnMatPC2PV,1));
            jp = ceil(rand(1)*size(ConnMatPC2PV,2));
            if ConnMatPC2PV(ip,jp) && ~ConnMatPV2PC(jp,ip)
                foundIt = 1;
            end
        end
        foundIt = 0;
        while ~foundIt
            ii = ceil(rand(1)*size(ConnMatPC2PV,1));
            ji = ceil(rand(1)*size(ConnMatPC2PV,2));
            if ~ConnMatPC2PV(ii,ji) && ConnMatPV2PC(ji,ii)
                foundIt = 1;
            end
        end
        % move the projection of the PC to form reciprocal with the PV:
        ConnMatPC2PV(ip,jp) = 0;
        ConnMatPC2PV(ii,ji) = 1;
        ctr = ctr - 1 ;
    end
    % Set indices for ease of mind:
    pc = size(PCsomata,1);
    pv = size(PCsomata,1) + size(PVsomata,1);
    cb = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1);
    cr = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1) + size(CRsomata,1);
    
    % Populate the final all-to-all connectivity matrix:
    AllConnMat = zeros(nAllCells);
        
    %         % Pyramidals to all types of interneurons:
    %         AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
    % PCs to PVs
    AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
    %     AllConnMat(1:pc,pv+1:end) = 1; % Connected to all interneurons: too much. maybe less connectivity?
    AllConnMat(1:pc,pv+1:end) = rand(nPC, nAllCells-pv) > 0.9; % Connected to all interneurons: too much. maybe less connectivity?
    % PVs connect to all other PVs + autapses (need to gap junctions?)
    AllConnMat(pc+1:pv,pc+1:pv) = 1;
    % PVs to PCs based on above connectivity (Yuste 2011):
    AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
    % CB only connect to PC (Xenia) na to psa3w...
    AllConnMat(pv+1:cb,1:pc) = ConnMatCB2PC;
    % CR only connect to PC and CB (Xenia) na to psa3w...
    AllConnMat(cb+1:cr,1:pc) = 1;
    AllConnMat(cb+1:cr,pv+1:cb) = 1;
    
    
    imagesc(AllConnMat)
    
    AllWeightMat = ones(size(AllConnMat));
    
    
    
    
    
%     connProbRange = -0.5:0.1:0.5;
    connProbRange = ones(1,11)*0.9;
    PC2PC_str = zeros(nPC,nPC,length(connProbRange));
    PC2PC_rnd = zeros(nPC,nPC,length(connProbRange));
    for i=1:length(connProbRange)
        factorOverallProb = 1+connProbRange(i);
        % Pyramidals connect to all
        
        connProbsPC2PC = connProbsPC2PC*factorOverallProb;
        recipProbsPC2PC = recipProbsPC2PC *factorOverallProb;
        NincomingProbs = NincomingProbs *factorOverallProb;
        NoutgoingProbs = NoutgoingProbs *factorOverallProb;
        
        MAXPC2PC = [];
        CCDIFF = [];
        CCMAX = 0;
        tic;
        for manytimes=1:200
            manytimes
            PC2PC_str(:,:,i) = create_graph_CN(distPC2PC,connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
                recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
            prob_conn_rnd_ind = sum(sum(PC2PC_str(:,:,i))) / (numel(PC2PC_str(:,:,i))-nPC); % joint probability
            PC2PC_rnd(:,:,i) = create_graph_WS(distPC2PC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
            [mcccs,~,cccs] = clust_coeff( PC2PC_str(:,:,i) );
            [mcccr,~,cccr] = clust_coeff( PC2PC_rnd(:,:,i) );
            CCDIFF(manytimes) = mcccs - mcccr;
            if(CCMAX < CCDIFF(manytimes))
                disp(sprintf('Replacing %f with %f',CCMAX,CCDIFF(manytimes)))
                CCMAX = CCDIFF(manytimes);
                MAXPC2PC(:,:,1) = PC2PC_str(:,:,i);
                MAXPC2PC(:,:,2) = PC2PC_rnd(:,:,i);
            end
        end
        runtime = tic-toc
        max(CCDIFF)
        PC2PC_str(:,:,i) = MAXPC2PC(:,:,1);
        PC2PC_rnd(:,:,i) = MAXPC2PC(:,:,2);
    end


    CCDIFF = [];
    for i=1:length(connProbRange)
        [mcccs,~,cccs] = clust_coeff( PC2PC_str(:,:,i) );
        [mcccr,~,cccr] = clust_coeff( PC2PC_rnd(:,:,i) );
        CCDIFF(i) = mcccs - mcccr;
    end
    
    for stateID = 1:size(PC2PC_str,3)
        PC2PC = [];
        PC2PC(:,:,1) = PC2PC_str(:,:,stateID);
        PC2PC(:,:,2) = PC2PC_rnd(:,:,stateID);

        % Find nearest neighbors of Pyramidals only
        [CNi_str,CNo_str] = m_commonNeighbors(PC2PC(:,:,1));
        [CNi_rnd,CNo_rnd] = m_commonNeighbors(PC2PC(:,:,2));
        mergedCN_str = [CNi_str + CNo_str] .* PC2PC(:,:,1);
        mergedCN_rnd = [CNi_rnd + CNo_rnd] .* PC2PC(:,:,2);
        mergedCN_strT = mergedCN_str';
        mergedCN_str(logical(tril(ones(size(mergedCN_str)),-1))) = mergedCN_strT(logical(tril(ones(size(mergedCN_str)),-1)));
        mergedCN_rndT = mergedCN_rnd';
        mergedCN_rnd(logical(tril(ones(size(mergedCN_rnd)),-1))) = mergedCN_strT(logical(tril(ones(size(mergedCN_rnd)),-1)));
        mergedCN_strNaN = mergedCN_str;
        mergedCN_strNaN(mergedCN_strNaN==0) = NaN;
        mergedCN_rndNaN = mergedCN_rnd;
        mergedCN_rndNaN(mergedCN_rndNaN==0) = NaN;
        
        %perform affinity propagation algorithm:
        Sid=1; % depricated variable from another clustering algo..
        
        % try to force different no of cluster to get as many clusters as you
        % can from the populations with high min cells per cluster:
        srtClNo = 10;
        setFlg = 1;
        while setFlg
            for i=1:20
                [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str,srtClNo);
                NC_str(Sid) = length(unique(idx_str));
                [~,~,labels_str(:,Sid)] = unique(idx_str);
                if(min(histc(labels_str(:,Sid),1:NC_str(Sid))') > 6)
                    setFlg = 0;
                    break;
                end
                
            end
            srtClNo = srtClNo - 1;
        end
        %     [idx_str,~,~,~]=apclustermex(mergedCN_str, nanmedian(mergedCN_strNaN,2)*medianFactor_str );
        cellsPerCluster_str = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
        figure('name','Structured')
        plot(cellsPerCluster_str{:})
        
        % try to force different no of cluster to get as many clusters as you
        % can from the populations with high min cells per cluster:
        srtClNo = 10;
        setFlg = 1;
        while setFlg
            for i=1:20
                [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd,srtClNo);
                NC_rnd(Sid) = length(unique(idx_rnd));
                [~,~,labels_rnd(:,Sid)] = unique(idx_rnd);
                if(min(histc(labels_rnd(:,Sid),1:NC_rnd(Sid))') > 6)
                    setFlg = 0;
                    break;
                end
                
            end
            srtClNo = srtClNo - 1;
        end
        
        %how many cells in each cluster?
        cellsPerCluster_rnd = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
        figure('name','Random')
        plot(cellsPerCluster_rnd{:})
        
        % run for the structured network and for the random network as
        % well:
        
        RUNS_str = {};
        RUNS_rnd = {};
        
        %     export network parameters for that experiment:
        cd(mypath);
        system(sprintf('mkdir experiment_%d',aa));
        exportNetworkPositions(PCsomata,PVsomata,sprintf('%sexperiment_%d/',mypath,aa));
        % Export pyramidal's clusters ID:
        exportNetworkCluster([idx_str,idx_rnd],sprintf('%sexperiment_%d/',mypath,aa));
        
        save(sprintf('%sexperiment_%d/experiment.mat',mypath,aa),'-v7.3');

        cellsPerCluster_str(stateID) = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
        cellsPerCluster_rnd(stateID) = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
        Clusters_str(1:nPC,stateID) = labels_str(:,Sid);

        stimCellsPerCluster = min(min(cellsPerCluster_str{1,stateID}),min(cellsPerCluster_rnd{1,stateID}));
        
    
        for t=1:max(Clusters_str(1:nPC,stateID))
            fprintf('@experiment %d, cluster %d... ',aa, t);
            % Isolate one cluster:
            [r,~,~] = find(Clusters_str(1:nPC,stateID) == t);
            StimVect_str = zeros(nPC,1);
            tmpr = randperm(length(r));
            StimVect_str(r(tmpr(1:stimCellsPerCluster))) = 1;  % 7==NUmber of total cells stimulated in each cluster.
            %         StimVect_str = ones(nPC,1); %stimulate WHOLE NETWORK
            Clusters_stim_str(1:nPC,stateID) = StimVect_str;
            % Autapses:
            AllConnMat_str = AllConnMat;
            AllConnMat_str(1:pc,1:pc) = PC2PC(:,:,1);
            AllConnMat_str(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
            %         AllConnMat_str(1:pc,1:pc) = 1;
            
            AllWeightMat_str = AllWeightMat;
            
            % Kayaguchi : more potent synapses in reciprocal pairs
            % Perin: more potent synapses with more clustering coeff:
            % We try to combine them
            [~,~,cccs] = clust_coeff( AllConnMat_str(1:pc,1:pc) );
            cccs = cccs-min(cccs)
            cccs = cccs/max(cccs)
            weightsexp = exp( [0:0.001:0.3]);
            weightsexp = max(weightsexp) - weightsexp
            weightsexp = weightsexp(end:-1:1);
            plot(weightsexp)
            
            for c = 1:nPC
%                 AllWeightMat_str(1:nPC,c) = 1 + (cccs(c)*1.5);
                    AllWeightMat_str(1:nPC,c) = 1 + (cccs(c)*weightsexp(c)*1.5);
            end
            
            %         AllWeightMat_str = AllWeightMat_str .* AllConnMat_str ;
%             AllWeightMat_str = ones(size(AllConnMat_str));
            
            %  ---- Export parameters to NEURON ----
            exportNetworkParameters(AllCells,AllConnMat_str,AllWeightMat_str,sprintf('%sexperiment_%d/',mypath,aa));
            exportStimulationParameters(AllCells,StimVect_str,sprintf('%sexperiment_%d/',mypath,aa));
            
            fprintf('Running NEURON... ');
            unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
            tic;
            VCLAMP = 0;
            EXPERIMENT =1;  % 0 = random, 1= structured gia to $#$#% NEURON
            SIMPLIFIED = 1;
            PARALLEL = 1;
            if (SIMPLIFIED)
                neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
                    '%smechanism_simple/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "EXPID=%d" '...
                    '-c "AA=%d" '...
                    '-c "VCLAMP=%f" '...
                    'final.hoc'], ...
                    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,stateID,aa, VCLAMP);
            else
                neuroncmd = sprintf(['mpirun -n 24 -mca btl ^openib '...
                    '%smechanism_complex/x86_64/special -mpi -nobanner '...
                    '-c "PARALLEL=%d" '...
                    '-c "SIMPLIFIED=%d" '...
                    '-c "CLUSTER_ID=%d" '...
                    '-c "EXPERIMENT=%d" '...
                    '-c "EXPID=%d" '...
                    '-c "AA=%d" '...
                    '-c "VCLAMP=%f" '...
                    'final.hoc'], ...
                    mypath2,PARALLEL,SIMPLIFIED,t,EXPERIMENT,stateID,aa, VCLAMP);
            end
            [nrnStatus,nrnCmdOut{t}] = unix(neuroncmd);
            runtime = toc
            
            if( ~findstr('Success',nrnCmdOut{t}))
                continue;
            end
            
            fprintf('DONE!\n');
        end % for t
    end
%% Load cluster info:
mypath = '/srv/userdata/HOMES/stefanos/Desktop/prefrontal-micro/experiment/network/';
mypath2 = '~/Desktop/prefrontal-micro/';
load(sprintf('%sexperiment_%d/experiment.mat',mypath,aa));
idx_str = load(sprintf('%sexperiment_%d/networkPyramidalClustersStructured.txt',mypath,aa));
idx_rnd = load(sprintf('%sexperiment_%d/networkPyramidalClustersRandom.txt',mypath,aa));
[~,~,labels_str] = unique(idx_str);
[~,~,labels_rnd] = unique(idx_rnd);
Sid = 1;
exprun = 1;
pc = 75;
stimend = 1500;
NC_str = length(unique(labels_str));
NC_rnd = length(unique(labels_rnd));
PCsomata = load(sprintf('%sexperiment_%d/networkPositionsPyramidals.txt',mypath,aa));

%%
fprintf('Loading runs...');
tic;
%choose structured or random:
% expString = 'RND';
expString = sprintf('STR_%d',stateID);

PCcells_str = {};
PCcells_str_spk = {};
PCcells_str_i = {};
PCcells_str_dist= {};
for t=1:max(Clusters_str(1:nPC,stateID))
    for ru = exprun
        for c=1:pc
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            mycell.position = PCsomata(c,1:3);
%             mycell.clusterID = Clusters_str(c,connProbSelect);
            PCcells_str{c,ru}=mycell;
%             PCcells_str_spk{c,ru} = load(sprintf('%sexperiment_%d/%s/%d_%d_%d_st.txt',mypath,aa,expString,t,c-1,ru-1));
%             PCcells_str_dist{c,ru} = load(sprintf('%sexperiment_%d/%s/%d_%d_%d_DIST.txt',mypath,aa,expString,t,c-1,ru-1));
        end
        
        for c=pc+1:pv
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            PVcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
        for c=pv+1:cb
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            CBcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
        for c=cb+1:cr
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
           CRcells_str{c,ru}=mycell;%.hasPersistent(1000,8,mycell.tstop-1001); % paper?
        end
    end
    RUNS_str{1,t} = PCcells_str(:,:);
end

runtime = toc
fprintf('DONE!\n');

%% Plot
% close all;

ru = 1 ;% what run of this cluster

FFS = [];
FFPA = [];
for t=1:max(Clusters_str(1:nPC,stateID))
    for c=1:nPC
%         comprobj(c) = RUNS_str{1,t}{c,ru};
        FFS(t,c) = RUNS_str{1,t}{c,ru}.spike_count(500,1500) ;
        FFPA(t,c) = RUNS_str{1,t}{c,ru}.spike_count(1500,tstop) / 3.5;
    end
end
% figure();
% plot(FFS');hold on;
figure();
plot(FFPA');hold on;
plot(std(FFPA),'k','linewidth',2);hold on;

% tmp=find(StimVect_str);
% for c = 1:length(tmp)
% plot(tmp(c),RUNS_str{1,t}{tmp(c),ru}.spike_count(1500,tstop) / 3.5,'+r');hold on;
% end
%% Firing Frequency:
close all;
PCFF_mean = [];

for dataset=1
    % load dataset:
    if dataset == 1
        NC = max(Clusters_str(1:nPC,stateID));
        RUNS = RUNS_str;
    else 
        NC = NC_rnd(Sid);
%         NC = NC_str(Sid);
        RUNS = RUNS_rnd;
    end



%     % if in random network, replace the clusterID with that of the
%     % structured, becuse we dont have real clusters in random condition
%     if dataset == 2
%        for clu = 1:NC_str
%            for ru = 1:exprun % what run of this cluster
%                 for c=1:pc
%                     RUNS{1,clu}{c,ru}.clusterID = RUNS_str{1,clu}{c,ru}.clusterID;
%                 end
%            end
%        end
%     end
    
    
    
    FF = [];
    TMP=[];
    for ru = exprun % what run of this cluster
        FFPerClust=cell(NC(Sid));
        PAPerClust=zeros(NC(Sid));

        for stimclu =1:NC(Sid)
            for cl=1:NC(Sid)
                for c=1:pc
                    if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
%                         if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                            %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                            responce = RUNS{1,stimclu}{c,ru}.mv((stimend*10)+1:end);
                            FFPerClust{cl,stimclu} = [FFPerClust{cl,stimclu} ,length(find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) )) / ((RUNS{1,stimclu}{c,ru}.tstop-(stimend+1))/1000) ];
%                         end
                    end
                end
            end
        end
        %     PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
        %     FFPC(:,:,ru) = cellfun(@(x,i) i*(nanmean(x)),FFPerClust,mat2cell(cellfun(@(x) ~isempty(x), FFPerClust),ones(NC(Sid),1),ones(NC(Sid),1)) );


        TMP = cellfun(@(x) nanmean(x),FFPerClust );
        TMP(~isfinite(TMP))  = 0;
%         FF(:,:,ru) = TMP./ ((RUNS{1,stimclu}{c,ru}.tstop-(stimend+1))/1000)  ; % Get Freq, no num of spikes
        FF(:,:,ru) = TMP;
    end

    FFnonstim = FF;

    FFstim = [];
    for i=exprun
        FFstim = [FFstim;diag(FFnonstim(:,:,i)) ];
        FFnonstim(:,:,i) = FFnonstim(:,:,i).*(~eye(NC(Sid)));
    end
    FFnonstim =FFnonstim(find(FFnonstim))  ;

    % non stimulated
    mean(FFnonstim)
    std(FFnonstim)

    % in stimulated cluster:
    mean(FFstim)
    std(FFstim)

    % Per Cluster:
    if dataset == 1
        PCFF_mean(:,:,1) = mean(FF,3);
    else 
        PCFF_mean(:,:,2) = mean(FF,3);
    end


    % ttest between the same cluster:
    for i=1:NC(Sid)
        for j=1:NC(Sid)
            [h(i,j),p(i,j)] = ttest2(reshape(FF(i,:,:),1,[]),reshape(FF(j,:,:),1,[]))
        end
    end

end% for datasets

ccm = jet( ceil(max(max(max(PCFF_mean)))*100)  ) ; 
% ccm = ccm(1:100:end,:);
figure('Name','Structured Firing Frequency' );
imagesc(PCFF_mean(:,:,1));
colormap(ccm(1:round(max(max(PCFF_mean(:,:,1)))*100),:));
colorbar;
set(gcf,'Color',[1,1,1]);

% figure('Name','Random Firing Frequency' );
% imagesc(PCFF_mean(:,:,2));
% colormap(ccm(1:round(max(max(PCFF_mean(:,:,2)))*100),:));
% colorbar;
% set(gcf,'Color',[1,1,1]);
% 
% figure('Name','Histo Firing Frequency' );
% ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,1),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
% ch = get(ph,'Children');
% set(ch,'FaceAlpha',0.6, 'FaceColor', 'b') ;
% ph = bar( 0:0.5:ceil(max(max(max(PCFF_mean)))), histc(reshape(PCFF_mean(:,:,2),1,[]),0:0.5:ceil(max(max(max(PCFF_mean)))) )' );hold on;
% ch = get(ph,'Children');
% set(ch,'FaceAlpha',0.6, 'FaceColor', 'r') ;
% set(gcf,'Color',[1,1,1]);
% 
% figure;
% bar(cpc);
% set(gcf,'Color',[1,1,1]);
%% load morphologies:
% addpath('~/Documents/MATLAB/Dendrites/');
% addpath(genpath('~/Documents/MATLAB/Dendrites/TREES1.15'));
% morphPath = '/srv/userdata/HOMES/stefanos/Documents/MATLAB/Dendrites/' ;
% 
% % Load morphologies names from Smith's lab
% morphDir = 'SmithMorphologies';
% names = dir([morphPath,morphDir,'/']);
% names = {names(~[names.isdir]).name};
% names = names(cellfun(@(x) strcmp(x(end-3:end),'.swc'), names));
% 
% % Load all morphologies as ncell objects:
% for i=1:length(names)
%     i
%     PC(i) =  ncell([],1,[morphPath,morphDir,'/',names{i}]);
% end


%% Plot
% close all;

ru = 1 ;% what run of this cluster
Incommers = {};
for c=1:nPC
    Incommers{c} = find(AllConnMat_str(1:nPC,c));
end
tmpArr=[];
tmpMax=[];
tmpEachArr={};
for c=1:nPC
    tmpEachArr{c}=[];
    for i=1:length(Incommers{c})
        tmpArr = [tmpArr,PCcells_str_dist{Incommers{c}(i),ru}];
        tmpEachArr{c} = [tmpEachArr{c},PCcells_str_dist{Incommers{c}(i),ru}];
    end
    tmpMax(c) = max(tmpEachArr{c});
end
% tmpMin = min(tmpArr);
% tmpMax = max(tmpArr);
cm = jet( ceil(max(tmpArr)) - floor(min(tmpArr)) );
cm(1,:) = [1,1,1];
gm = gray(200);
gm = gm(end:-1:100,:); % overhead; reversed: not starting from black
colmap = [cm;gm];

figure();
for c=1:nPC
    colormap(colmap);
    c
% %     subplot('Position',[0.05,0.05,0.9,0.1]);
% subplot('Position',[0.05,0.55,0.9,0.4]);
% %     The raster plot below is VERY heavy...
% %     for syn = 1:length(Incommers{c})
% %         for l=1:size(PCcells_str_spk{Incommers{c}(syn),ru},1)
% %             plot([PCcells_str_spk{Incommers{c}(syn),ru}(l),PCcells_str_spk{Incommers{c}(syn),ru}(l)],[syn-1,syn],'k');hold on;
% %         end 
% %         HIST(syn) = PCcells_str_dist{Incommers{c}(syn),ru};
% %     end
%     imagesc( repmat( sum( incommingSynEvents{c} ) , length(Incommers{c}), [] ) + ( ceil(max(tmpArr)) - floor(min(tmpArr)) ) );
%     axis([0,tstop,0,length(Incommers{c})])
%     set(gca,'ytick',[]);
%     set(gca,'xtick',[]);
%     hold on;
    
    subplot('Position',[0.05,0.55,0.9,0.4]);
    plot(RUNS_str{1,t}{c,ru}.mv);
    axis([0,tstop*10,-90,60]);
    
    subplot('Position',[0.05,0.2,0.9,0.3]);
    incommingSynEvents{c}=zeros(length(Incommers{c}),tstop);
    
    for syn = 1:length(Incommers{c})
        for l=1:size(PCcells_str_spk{Incommers{c}(syn),ru},1)
            incommingSynEvents{c}(Incommers{c}(syn),round(PCcells_str_spk{Incommers{c}(syn),ru}(l))) = 1;
        end 
    end
    axis([0,tstop,0,length(Incommers{c})])
    ylabel('Incomming Synaptic Events');
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    imagesc(incommingSynEvents{c});
%     colormap( cm(1:round(tmpMax(c)) - floor(min(tmpArr)),:) );
    disp(sprintf('Mean distance from soma is: %f\n' , mean(tmpEachArr{c}) ))
    disp(sprintf('Incomming connections are: %f\n' , length(Incommers{c}) ))
    
    pause;
    cla;
end


%% detect UP states:

for c=1:nPC
    MUA(c,:) = RUNS_str{1,2}{c,1}.mv';
end
MUA = mean(MUA,1);

% get UP states more than 0.1 sec (for now):
[S,D,UPs]=findUPstates(MUA, 2, 3, RESTING, 1000 ) ;

% % Truncate all by the shorter UP state:
% UP_st = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,UPs)))',UPs,'UniformOutput',false ) ) ;
% 
% UP_sa = mean(UP_st,1) ;
% 
% plot(UP_st', 'color', [0.9,0.9,0.9]);hold on;
% plot(UP_sa', 'k');


cl = lines(7);
IVu=[];
figure();
for t=2:8
    for i=1:length(D)
        VCi{i,t} = PCcells_str_i{t-1,1}.mv(S(i):S(i)+D(i));
    end
    % Truncate all by the shorter UP state:
    VCi_avr = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,VCi(:,t))))',VCi(:,t),'UniformOutput',false ) ) ;
    VCi_avr = mean(VCi_avr,1) ;
    
    plot(VCi_avr, 'color',cl(t-1,:));hold on;
    plot([0,length(VCi_avr)],[min(VCi_avr),min(VCi_avr)], 'color',(cl(t-1,:)+0.5)/norm(cl(t-1,:)+0.5),'linewidth',2);hold on;
    lstr{t-1} = sprintf('Vclamp @%.1f',RESTING + test_vclamp(t-1));
    IVu(t-1,:) = min(VCi_avr) ./ RESTING + test_vclamp(t-1);
end
legend(lstr);




MUA_DOWN=[];
IVd=[];
figure();
for t=2:8
    MUA_DOWN(t-1,:) = PCcells_str_i{t-1,1}.mv(21350:24780);
    plot(MUA_DOWN(t-1,:), 'color',cl(t-1,:));hold on;
    plot([0,length(MUA_DOWN(t-1,:))],[min(MUA_DOWN(t-1,:)),min(MUA_DOWN(t-1,:))], 'color',(cl(t-1,:)+0.5)/norm(cl(t-1,:)+0.5),'linewidth',2);hold on;
    IVd(t-1,:) = min(MUA_DOWN(t-1,:)) / RESTING + test_vclamp(t-1);
end
legend(lstr);


figure();
plot(IVu);hold on;
max(diff(IVu))
plot(IVd,'r');hold on;
max(diff(IVd))
%%

plot( PCcells_str{c,ru}.mv & ((PCcells_str{c,ru}.mv-RESTING)>5) );
hold on;
plot(PCcells_str{c,ru}.mv, 'r')

%%

t=1; % Initial without vclamps
VCLAMP = 0;
GPYID = 0;

fprintf('Running NEURON... ');
unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
EXPERIMENT = 1;  % 0 = random, 1= structured gia to $#$#% NEURON
%         [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpiexec -np 8 %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
[nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" -c "VCLAMP=%f"  -c "GPYID=%d" final.hoc',mypath2,t,EXPERIMENT,aa, VCLAMP, GPYID));
runtime = toc



%% detect UP states; Current clamp:

for c=2:nPC
    MUA(c,:) = RUNS_str{1,1}{c,1}.mv';
end
MUA = mean(MUA,1);

% get UP states more than 0.1 sec (for now):
[S,D,UPs]=findUPstates(MUA, 2, 3, RESTING, 1000 ) ;


for i=1:length(D)
    VC{i} = RUNS_str{1,t}{1,1}.mv(S(i):S(i)+D(i));
end
% Truncate all by the shorter UP state:
VC_avr = cell2mat( cellfun(@(x) x(1:min(cellfun(@length,VC)))',VC,'UniformOutput',false )' ) ;
VC_avr = mean(VC_avr,1) ;

plot(VC_avr)

Rin_active = (VC_avr(1) - min(VC_avr)) / (-0.2) ;

