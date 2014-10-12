
for aa = 1:1
    close all; clearvars -except aa; clc;
    aa
    addpath(genpath('~/Documents/MATLAB/'));
    addpath(genpath('~/Documents/GitHub/prefrontal-micro/'));
    addpath(genpath('~/Documents/GitHub/neurocommitto/'));
    mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
    %cd(mypath);
    
    exprun=2;
    nPC = 75;%216;
    clstr = 1; % Cluster type
    
    connBinsPC2PC = 20:30:500;
    connProbsPC2PC = 0.25 .* exp(-0.006 * connBinsPC2PC);
    connBinsPV2PC = 0:20:500;
    connProbsPV2PC = linspace(1,0,length(connBinsPV2PC)) * 0.7;
    connBinsCB2PC = 0:20:500;
    connProbsCB2PC = normpdf(connBinsCB2PC,150,90) *10;
    
    recipBinsPC2PC = 20:30:500;
    recipProbsPC2PC = 0.12 .* exp(-0.006 * recipBinsPC2PC);
    % Update probabilities from Perin et al. (Figure S4)
    NconnBins = [0,1,2,3];
    NincomingProbs = [0.1, 0.2, 0.25, 0.4] ;
    NoutgoingProbs = [0.12, 0.25, 0.24, 0.2] ;
    
    %  ---- Initialize Pyramidal cells ----
    PCsomata = CreateRandomNetwork(nPC, 200, 3);
    distPC2PC = generateDistanceMat(PCsomata', 0);
    
    % Pyramidals connect to all
    %     glProbConn = 0.9;
    PC2PC = zeros(nPC,nPC,2);
    %     PC2PC(:,:,5) = create_graph_WS(distPC2PC,0.15,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    %     PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.16,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    %     PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.2,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    %     PC2PC(:,:,2) = create_graph_DD(distPC2PC,0.98,connBinsPC2PC, connProbsPC2PC);
    %     PC2PC(:,:,1) = create_graph_CN(distPC2PC,0.07, connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
    %         recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    

    PC2PC(:,:,1) = create_graph_CN(distPC2PC,connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    
    prob_conn_rnd_joint = sum(sum(tril(PC2PC(:,:,1) | PC2PC(:,:,1)') )) / ((numel(PC2PC(:,:,1))-nPC)/2); % individual probability
    prob_conn_rnd_int = sum(sum(PC2PC(:,:,1))) / (numel(PC2PC(:,:,1))-nPC); % joint probability
    
    PC2PC(:,:,2) = create_graph_WS(distPC2PC,prob_conn_rnd_joint,'joint'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    
    %     PC2PC(:,:,2) = create_graph_DD(distPC2PC,connBinsPC2PC, connProbsPC2PC);
    
    [mcccs,~,ccc] = clust_coeff( PC2PC(:,:,1) );
    [mcccr,~,ccc] = clust_coeff( PC2PC(:,:,2) );
    [mcccr,mcccs]
    
    
    %     tempPC2PC(:,:) = PC2PC(:,:,clstr);
    
    % you constrain the probability of connection to avoid
    % overconnectivity, but are the occuring networks retain their
    % distinctive clustering properties?
    
    % routines to calculate and plot the degrees and clustering
    % distribution:
    %     c_m = hot(6);
    %     for i=1:5
    %         p_h(i) = plot(1:5:100,histc(degrees(PC2PC(:,:,i)),1:5:100),'color',c_m(i,:),'linewidth',3);hold on;
    %
    %         set(get(get(p_h(i),'Annotation'),'LegendInformation'),...
    %     'IconDisplayStyle','on');
    %     end
    %     legend('Perin','Distance','Clustered','intermediate','Random');
    
    
    % ---- affinity propagation similarity measure: ----
    
    % Find nearest neighbors of Pyramidals only
    [CNi_str,CNo_str] = m_commonNeighbors(PC2PC(:,:,1));
    [CNi_rnd,CNo_rnd] = m_commonNeighbors(PC2PC(:,:,2));
    mergedCN_str = [CNi_str + CNo_str] .* PC2PC(:,:,1);
    mergedCN_rnd = [CNi_rnd + CNo_rnd] .* PC2PC(:,:,2);
    
    %perform affinity propagation algorithm:
    Sid=1;
    [idx_str,~,~,~]=apclustermex(mergedCN_str,repmat(median(mergedCN_str(find(mergedCN_str)))*0.5,[1,nPC]) );
    [idx_rnd,~,~,~]=apclustermex(mergedCN_rnd,repmat(median(mergedCN_rnd(find(mergedCN_rnd)))*0.5,[1,nPC]) );
    NC_str(Sid) = length(unique(idx_str));
    [~,~,labels_str(:,Sid)] = unique(idx_str);
    NC_rnd(Sid) = length(unique(idx_rnd));
    [~,~,labels_rnd(:,Sid)] = unique(idx_rnd);
    %     [NC,labels,Sid] = runAffinityPropagation(onlyCN);
    
    
    %how many cells in each cluster?
    cellsPerCluster_str = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
    cellsPerCluster_rnd = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
    figure('name','Structured')
    bar(cellsPerCluster_str{:})
    figure('name','Random')
    bar(cellsPerCluster_rnd{:})
    
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
    
    
    
    PVsomata = CreateCubeNetworkPV(250, nPV); % 226.6 per mm squared (?!?) paper?
    CBsomata = CreateCubeNetworkPV(200, nCB);
    CRsomata = CreateCubeNetworkPV(200, nCR);
    
    %     scatter3(PVsomata(:,1),PVsomata(:,2), PVsomata(:,3), 30, rand(1,3), 'fill', 'o');hold on;
    %     scatter3(CBsomata(:,1),CBsomata(:,2), CBsomata(:,3), 30, rand(1,3), 'fill', 'o');hold on;
    %     scatter3(CRsomata(:,1),CRsomata(:,2), CRsomata(:,3), 30, rand(1,3), 'fill', 'o');hold on;
    %     axis equal;
    
    
    % Distance of each interneuron from Pyramidal and connect based on
    % probability from Yuste 2011:
    distPV2PC = distancePV2PC(PVsomata,PCsomata);
    distCB2PC = distancePV2PC(CBsomata,PCsomata);
    
    ConnMatPV2PC = connectPV2PC(distPV2PC,connBinsPV2PC,connProbsPV2PC);
    ConnMatCB2PC = connectCB2PC(distCB2PC,connBinsCB2PC,connProbsCB2PC);
    
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
    % Set indies for ease of mind:
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
    AllConnMat(1:pc,pv+1:end) = 1; % Connected to all interneurons
    % PVs connect to all other PVs + autapses (need to gap junctions?)
    AllConnMat(pc+1:pv,pc+1:pv) = 1;
    % PVs to PCs based on above connectivity (Yuste 2011):
    AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
    % CB only connect to PC (Xenia) na to psa3w...
    AllConnMat(pv+1:cb,1:pc) = ConnMatCB2PC;
    % CR only connect to PC and CB (Xenia) na to psa3w...
    AllConnMat(cb+1:cr,1:pc) = 1;
    AllConnMat(cb+1:cr,pv+1:cb) = 1;
    
    AllWeightMat = AllConnMat;
    
    % run for the structured network and for the random network as
    % well:
    
    RUNS_str = {};
    RUNS_rnd = {};
    
    for t=1:NC_str(Sid)
        fprintf('@experiment %d... ',t);
        % Isolate one cluster:
        [r,~,~] = find(labels_str(:,Sid) == t);
        StimVect_str = zeros(nPC,1);
        %choose randomly three cells to stimulate:
        %         rp = randperm(length(r));
        %         StimVect(r(rp(1:3))) = 1;
        StimVect_str(r) = 1;
        %         if(cellsPerCluster{:}(t) > 3)
        % tempPC2PC(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1; % autapses in Pyramidals
        
        %         % Make all-to-all connections inside clusters!
        %         for ec = 1:NC(Sid)
        %             clid = find(labels(:,Sid) == ec);
        %             tempPC2PC(clid,clid) = 1;
        %         end
        AllConnMat_str = AllConnMat;
        AllConnMat_str(1:pc,1:pc) = PC2PC(:,:,1);
        AllConnMat_str(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
        
        AllWeightMat_str = AllWeightMat;
        
        % Kayaguchi MONO GIA TA RECIPROCAL!
        % OXI TA STO IDIO CLUSTER
        %         for ec = 1:NC(Sid)
        %             clid = find(labels(:,Sid) == ec);
        %             AllWeightMat(clid,clid) = 1.2;
        %         end

        %  ---- Export parameters to NEURON ----
        exportNetworkParameters(AllCells,AllConnMat_str,AllWeightMat_str,mypath);
        exportStimulationParameters(AllCells,StimVect_str,mypath);
        
        fprintf('Running NEURON... ');
        tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
        [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib ../../mechanism_complex/x86_64//special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" final.hoc',t));
        runtime = toc
        
        if( ~findstr('Success',nrnCmdOut{t}))
            continue;
        end
        
        
        fprintf('DONE!\n');
        
        %         else
        %             continue;
        %         end
    end % for t
    
    fprintf('Loading runs...');
    for t=1:NC_str(Sid)
        for ru = 1:exprun
            for c=1:pc
                mycell = ncell(load(sprintf('multi_core/%d_%d_%d.txt',t,c-1,ru-1)),10);
                mycell.clusterID = labels_str(c,Sid);
                mycell.position = PCsomata(c,1:3);
                PCcells{c,ru}=mycell.hasPersistent(1000,8,mycell.tstop-1001); % paper?
            end
%             for c=pc+1:pv
%                 mycell = ncell(load(sprintf('multi_core/%d_%d_%d.txt',t,c-1,ru-1)),10);
%                 mycell.position = PVsomata(c-pc,1:3);
%                 PVcells{c-pc,ru}=mycell;
%             end
%             for c=pv+1:cb
%                 mycell = ncell(load(sprintf('multi_core/%d_%d_%d.txt',t,c-1,ru-1)),10);
%                 mycell.position = CBsomata(c-pv,1:3);
%                 CBcells{c-pv,ru}=mycell;
%             end
%             for c=cb+1:cr
%                 mycell = ncell(load(sprintf('multi_core/%d_%d_%d.txt',t,c-1,ru-1)),10);
%                 mycell.position = CRsomata(c-cb,1:3);
%                 CRcells{c-cb,ru}=mycell;
%             end
        end
        RUNS_str{1,t} = PCcells(:,:);
%         RUNS_str{2,t} = PVcells(:,:);
%         RUNS_str{3,t} = CBcells(:,:);
%         RUNS_str{4,t} = CRcells(:,:);
%         Percentages(t) = sum(all(cellfun(@(x) x.persistentActivity,RUNS{1,t}))) * 100 / size(RUNS{1,t},2);
%         fprintf('Saving results... ');
%         save('batch_simplified_1_sameNoStimulated.mat','-v7.3');
    end
    fprintf('DONE!\n');
    
    save(sprintf('../../../%d_N%d_R%d_C%d.mat',aa,clstr,exprun,NC(Sid)),'-v7.3');
    % h5write('../../../CRUN.hdf5', '/dataset1', data);
    % hdf5read('test.hdf5', '/dataset1');
    
end
%%
for i=1:NC(Sid)
    clustID(:,i) = (labels(:,Sid) == i);
    CPA{i} = cellfun(@(x) x.persistentActivity,RUNS{1,i}(clustID(:,i),:));
    NCPA{i} = cellfun(@(x) x.persistentActivity,RUNS{1,i}(~clustID(:,i),:));
end

% generate color maps:
cm = rand(NC(Sid),3);
idcs= mat2cell(1:NC(Sid),[1],ones(NC(Sid),1));

figure();hold on;
cellfun(@(x,i) plot(histc(sum(x,1),0:1:nPC),'color',cm(i,:)), CPA,idcs,'uniformoutput',false);

%% Plot-a-cluster
close all;

clu = 1 ;% what cluster is stimulated
ru = 1;% what run of this cluster

CellsPerClust=zeros(NC(Sid),1);
PAPerClust=zeros(NC(Sid),1);


for cl=1:NC(Sid)
    figure(); hold on;
    set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
    for c=1:nPC
        if(RUNS{1,clu}{c,ru}.clusterID == cl)
            %                         figure;
            plot( RUNS{1,clu}{c,ru}.mv,'linewidth',1,'color',rand(1,3));
            CellsPerClust(cl) = CellsPerClust(cl)+1;
            if(RUNS{1,clu}{c,ru}.persistentActivity)
                PAPerClust(cl) = PAPerClust(cl) +1;
            end
        end
    end
    axis([0,length(RUNS{1,clu}{c,ru}.mv),-70,70 ]);hold on;
end

% bar graph
figure();
bar(CellsPerClust,'b');hold on;
bar(PAPerClust,'r');

% Extract some stats:

%% PA percentage
close all;

for ru = 1:10 % what run of this cluster
    
    CellsPerClust=zeros(NC(Sid));
    PAPerClust=zeros(NC(Sid));
    
    for clu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,clu}{c,ru}.clusterID == cl)
                    if(RUNS{1,clu}{c,ru}.persistentActivity)
                        PAPerClust(cl,clu) = PAPerClust(cl,clu) +1;
                    end
                end
            end
        end
    end
    PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
end

PAPSD = reshape(sum(PAP),[NC(Sid),10,1])';
mean(PAPSD)
std(PAPSD)

%% Firing Frequency:
close all;

for ru = 1:exprun % what run of this cluster
    FFPerClust=cell(NC(Sid));
    PAPerClust=zeros(NC(Sid));
    
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        responce = RUNS{1,stimclu}{c,ru}.mv(10001:end);
                        FFPerClust{cl,stimclu} = [FFPerClust{cl,stimclu} ,length(find( ([0;diff(sign(diff(responce)))<0;0] & [sign(responce)==1]) ))];
                    end
                end
            end
        end
    end
    %     PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
    %     FFPC(:,:,ru) = cellfun(@(x,i) i*(nanmean(x)),FFPerClust,mat2cell(cellfun(@(x) ~isempty(x), FFPerClust),ones(NC(Sid),1),ones(NC(Sid),1)) );
    
    
    TMP = cellfun(@(x) nanmean(x),FFPerClust );
    TMP(find(~isfinite(TMP)))  = 0;
    FFPC(:,:,ru) = TMP%./4  ; % why?
end

FFPCtmp = FFPC;

ffstimclu = [];
for i=1:exprun
    ffstimclu = [ffstimclu;diag(FFPCtmp(:,:,i)) ];
    FFPCtmp(:,:,i) = FFPCtmp(:,:,i).*(~eye(NC(Sid)));
end
FFPCtmp =FFPCtmp(find(FFPCtmp))  ;

% non stimulated
mean(FFPCtmp)
std(FFPCtmp)

% in stimulated cluster:
mean(ffstimclu)
std(ffstimclu)

%% Coefficient of Variation:
close all;

for ru = 1:exprun % what run of this cluster
    ISIPerClust=cell(NC(Sid));
    PAPerClust=zeros(NC(Sid));
    
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        ISIPerClust{cl,stimclu} = [ISIPerClust{cl,stimclu} ,diff(RUNS{1,stimclu}{c,ru}.spikes)'];
                    end
                end
            end
        end
    end
    %     PAP(:,:,ru) = (PAPerClust.*(~eye(length(PAPerClust))));
    %     FFPC(:,:,ru) = cellfun(@(x,i) i*(nanmean(x)),FFPerClust,mat2cell(cellfun(@(x) ~isempty(x), FFPerClust),ones(NC(Sid),1),ones(NC(Sid),1)) );
    
    
    %     TMP = cellfun(@(x) nanmean(x),ISIPerClust );
    %     TMP(find(~isfinite(TMP)))  = 0;
    ISIPC(:,:,ru) = ISIPerClust%./4  ;
end

% for non stimulated condition (per cluster):
ISIPCtmp = ISIPC;

for k = 1:exprun
    for i=1:NC(Sid)
        ISIPCtmp{i,i,k} = [];
    end
end

ISIinclust={};
for i=1:NC(Sid)
    ISIinclust{i}=[];
    for k=1:exprun
        ISIinclust{i} = [ISIinclust{i},ISIPCtmp{i,:,k}];
    end
end

bar(cellfun(@std,ISIinclust)./cellfun(@mean,ISIinclust))

% for the stimulated condition: (per cluster)
ISIPCdiag = {};
for i=1:NC(Sid)
    ISIPCdiag{i} = [];
    for k = 1:exprun
        ISIPCdiag{i} = [ISIPCdiag{i},ISIPC{i,i,k}] ;
    end
end

bar(cellfun(@std,ISIPCdiag)./cellfun(@mean,ISIPCdiag))

%% Cross correlation:
close all;
PerClust = cell(NC(Sid),NC(Sid),exprun,nPC);
for ru = 1:exprun % what run of this cluster
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
                        %                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        tmpSpikes = RUNS{1,stimclu}{c,ru}.spikes;
                        PerClust{cl,stimclu,ru,c} = tmpSpikes(tmpSpikes > 1001);
                    end
                end
            end
        end
    end
end

% Init SPIKY:
para.tmin = 0;
para.tmax = 5000;
para.dts =  10;
para.dtm = 10;
para.select_measures=[0 1     0 0     0 0     0 0];            % Select order of measures
% para.select_measures=[1];            % Select order of measures


% Check synchronicity of cells belonging to the stimulated cluster:
synchInsideClust={};
synchInsideClustclean={};
overal_dissim_stimulated={};
for i=1:NC(Sid)
    synchInsideClust{i}=[];
    for k=1:exprun
        synchInsideClust(i,k) = {reshape(PerClust(i,i,k,reshape(~cellfun(@isempty,PerClust(i,i,k,:)),1,[])),1,[])};
    end
    %concatinate all the cells of this cluster/run to one pool:
    % get idx of runs containing more than 1 cell firing after stimuly
    % (PA):
    goodIdx = cellfun(@(x) x>1,cellfun(@(x) length(x),synchInsideClust(i,:),'uniformoutput',false));
    synchInsideClustclean(i,1:sum(goodIdx)) = synchInsideClust(i,goodIdx) ;
    
    % feed SPIKY:
    for k=1:sum(goodIdx)
        results = SPIKY_no_plot_f_distances_MEX(synchInsideClustclean{i,k},para);
        overal_dissim_stimulated{i,k} = results.overall_dissimilarities ;
    end
    
    %     synchInsideClust{i,~cellfun(@isempty,synchInsideClust(i,:))}
end
% mean synchronicity of all the stimulated clusters:
mean_synch_stim = mean(cell2mat(overal_dissim_stimulated(cellfun(@(x) ~isempty(x), overal_dissim_stimulated))));
std_synch_stim = std(cell2mat(overal_dissim_stimulated(cellfun(@(x) ~isempty(x), overal_dissim_stimulated))));

% Check synchronicity of cells belonging to the non-stimulated cluster:
synchOutsideClust={};
syncnOutsideClustclean={};
overal_dissim_nostim={};
for i=1:NC(Sid)
    synchOutsideClust{i}=[];
    ctrk=1;
    for j=1:NC(Sid)
        for k=1:exprun
            if i ~= j
                synchOutsideClust(i,ctrk) = {reshape(PerClust(i,j,k,reshape(~cellfun(@isempty,PerClust(i,j,k,:)),1,[])),1,[])};
                ctrk = ctrk+1;
            end
        end
    end
    %concatinate all the cells of this cluster/run to one pool:
    % get idx of runs containing more than 1 cell firing after stimuly
    % (PA):
    goodIdx = cellfun(@(x) x>1,cellfun(@(x) length(x),synchOutsideClust(i,:),'uniformoutput',false));
    syncnOutsideClustclean(i,1:sum(goodIdx)) = synchOutsideClust(i,goodIdx) ;
    
    % feed SPIKY:
    for k=1:sum(goodIdx)
        results = SPIKY_no_plot_f_distances_MEX(syncnOutsideClustclean{i,k},para);
        overal_dissim_nostim{i,k} = results.overall_dissimilarities ;
    end
    
    %     synchInsideClust{i,~cellfun(@isempty,synchInsideClust(i,:))}
end
% mean synchronicity of all the stimulated clusters:
mean_synch_nostim = mean(cell2mat(overal_dissim_nostim(cellfun(@(x) ~isempty(x), overal_dissim_nostim))));
std_synch_nostim = std(cell2mat(overal_dissim_nostim(cellfun(@(x) ~isempty(x), overal_dissim_nostim))));


% Check synchronicity between stimulated and non stimulated clusters
% (pairwise):
% trains_a = {};
%     trains_b ={};
% for i=1:11
%     trains_a{i} = ceil(rand(ceil(rand(1)*20),1)*100) ;
%     trains_b{i} = ceil(rand(ceil(rand(1)*20),1)*100) ;
% end
% trains_a = trains_a';
% trains_b = trains_b';
% output = MNvanRossum(trains_a,trains_b,0.9,0.1)

syncnPairInClust={};
syncnPairInClustClean={};
syncnPairOutClust={};
syncnPairOutClustClean={};
overal_dissim_pairwise={};
for i=1:NC(Sid)
    for j=1:NC(Sid)
        % get spikes from stimulated cluster:
        if(i ~= j)
            %         fprintf('taking from ii %d\n',i)
            syncnPairInClust{j}=[];
            ctrk=1;
            for k=1:exprun
                syncnPairInClust(j,ctrk) = {reshape(PerClust(i,i,k,reshape(~cellfun(@isempty,PerClust(i,i,k,:)),1,[])),1,[])} ;
                ctrk = ctrk+1;
            end
            goodIdxIn = cellfun(@(x) x>1,cellfun(@(x) length(x),syncnPairInClust(j,:),'uniformoutput',false));
            syncnPairInClustClean(i,1:sum(goodIdxIn)) = syncnPairInClust(j,goodIdxIn) ;
            
            
            
            % get spikes from all other clusters in different cells:
            syncnPairOutClust{j}=[];
            ctrk=1;
            for k=1:exprun
                if i ~= j
                    syncnPairOutClust(j,ctrk) = {reshape(PerClust(i,j,k,reshape(~cellfun(@isempty,PerClust(i,j,k,:)),1,[])),1,[])};
                    ctrk = ctrk+1;
                end
            end
            
            goodIdxOut = cellfun(@(x) x>1,cellfun(@(x) length(x),syncnPairOutClust(j,:),'uniformoutput',false));
            syncnPairOutClustClean(i,1:sum(goodIdxOut)) = syncnPairOutClust(j,goodIdxOut) ;
            
            min_elements = [];
            min_rndperm_els=[];
            if( ~isempty(syncnPairOutClustClean) && ~isempty(syncnPairInClustClean) )
                min_elements = min(length(syncnPairInClustClean{1}), length(syncnPairOutClustClean{1}));
                tmpels = randperm(length(syncnPairInClustClean{1}));
                min_rndperm_els(1,:) = tmpels(1:min_elements);
                tmpels = randperm(length(syncnPairOutClustClean{1}));
                min_rndperm_els(2,:) = tmpels(1:min_elements);
                
                min_elements = min(length(syncnPairInClustClean{1}), length(syncnPairOutClustClean{1}));
                randperm(min_elements);
                overal_dissim_pairwise{i,j} = MNvanRossum(syncnPairInClustClean{1}(min_rndperm_els(1,:))',syncnPairOutClustClean{1}(min_rndperm_els(2,:) )',0.9,0.1)  ;
                %             fprintf('comparing with jj %d\n',j)
            end
            
        end
    end
end

overal_dissim_pairwise(cellfun(@isempty,overal_dissim_pairwise)) = {0};
overal_dissim_pairwise = cell2mat(overal_dissim_pairwise);

%%
scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 20, [0,0,1], 'fill', 'o');
axis equal; hold on;