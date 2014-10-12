
for aa = 1:1
    close all; clearvars -except aa; clc;
    aa
    addpath(genpath('~/Documents/MATLAB/'));
    addpath(genpath('~/Documents/GitHub/prefrontal-micro/'));
    addpath(genpath('~/Documents/GitHub/neurocommitto/'));
    mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';
    %cd(mypath);
    
    exprun=1;
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
    PC2PC = zeros(nPC,nPC,4);
%     PC2PC(:,:,5) = create_graph_WS(distPC2PC,0.15,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
%     PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.16,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
%     PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.2,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
%     PC2PC(:,:,2) = create_graph_DD(distPC2PC,0.98,connBinsPC2PC, connProbsPC2PC);
%     PC2PC(:,:,1) = create_graph_CN(distPC2PC,0.07, connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA 
%         recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    PC2PC(:,:,5) = create_graph_WS(distPC2PC,0.2,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.2,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.2,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    PC2PC(:,:,2) = create_graph_DD(distPC2PC,0.2,connBinsPC2PC, connProbsPC2PC);
    PC2PC(:,:,1) = create_graph_CN(distPC2PC,0.2, connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA 
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    
    tempPC2PC(:,:) = PC2PC(:,:,clstr);
    
    % you constrain the probability of connection to avoid
    % overconnectivity, but are the occuring networks retain their
    % distinctive clustering properties?
    
    % routines to calculate and plot the degrees and clustering
    % distribution:
    c_m = hot(6);
    for i=1:5
        p_h(i) = plot(1:5:100,histc(degrees(PC2PC(:,:,i)),1:5:100),'color',c_m(i,:),'linewidth',3);hold on;
        
        set(get(get(p_h(i),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','on'); 
    end
    legend('Perin','Distance','Clustered','intermediate','Random');
    
    
    % ---- affinity propagation similarity measure: ----
    
    % Find nearest neighbors of Pyramidals only
    [CNi,CNo] = commonNeighbors(tempPC2PC);
    onlyCN = [CNi + CNo] .* tempPC2PC;
    
    %perform affinity propagation algorithm:
    [NC,labels,Sid] = runAffinityPropagation(onlyCN);
    
    %how many cells in each cluster?
    cellsPerCluster = {histc(labels(:,Sid),1:NC(Sid))'};
    %     bar(cellsPerCluster{:})
    
    %
    for t=1:NC(Sid)
        fprintf('@experiment %d... ',t);
        % Isolate one cluster:
        [r,~,~] = find(labels(:,Sid) == t);
        StimVect = zeros(nPC,1);
        StimVect(r) = 1;
        
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
        
        
        % Populate the final all-to-all connectivity matrix:
        AllConnMat = zeros(nAllCells);
        % Set indies for ease of mind:
        pc = size(PCsomata,1);
        pv = size(PCsomata,1) + size(PVsomata,1);
        cb = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1);
        cr = size(PCsomata,1) + size(PVsomata,1) + size(CBsomata,1) + size(CRsomata,1);
        
        % tempPC2PC(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1; % autapses in Pyramidals
        tempPC2PC(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;
        
%         % Make all-to-all connections inside clusters!
%         for ec = 1:NC(Sid)
%             clid = find(labels(:,Sid) == ec);
%             tempPC2PC(clid,clid) = 1;
%         end
        
        AllConnMat(1:pc,1:pc) = tempPC2PC;
        
        
%         % Pyramidals to all types of interneurons:
%         AllConnMat(1:pc,pc+1:end) = 1; % Connected to all
        % PCs to PVs
        AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
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
        % Kayaguchi MONO GIA TA RECIPROCAL!
        % OXI TA STO IDIO CLUSTER
        for ec = 1:NC(Sid)
            clid = find(labels(:,Sid) == ec);
            AllWeightMat(clid,clid) = 1.2;
        end
        
        
        
        %  ---- Export parameters to NEURON ----
        exportNetworkParameters(AllCells,AllConnMat,AllWeightMat,mypath);
        exportStimulationParameters(AllCells,StimVect,mypath);
        
        fprintf('Running NEURON... ');
        [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
        if( ~findstr('Success',nrnCmdOut{t}))
            continue;
        end
        
        fprintf('Loading runs... ');
        for ru = 1:exprun
            for c=1:pc
                mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,ru-1)),10);
                mycell.clusterID = labels(c,Sid);
                mycell.position = PCsomata(c,1:3);
                PCcells{c,ru}=mycell.hasPersistent(1000,8,3999); % paper?
            end
            for c=pc+1:pv
                mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,ru-1)),10);
                mycell.position = PVsomata(c-pc,1:3);
                PVcells{c-pc,ru}=mycell;
            end
            for c=pv+1:cb
                mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,ru-1)),10);
                mycell.position = CBsomata(c-pv,1:3);
                CBcells{c-pv,ru}=mycell;
            end
            for c=cb+1:cr
                mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,ru-1)),10);
                mycell.position = CRsomata(c-cb,1:3);
                CRcells{c-cb,ru}=mycell;
            end
        end
        RUNS{1,t} = PCcells(:,:);
        RUNS{2,t} = PVcells(:,:);
        RUNS{3,t} = CBcells(:,:);
        RUNS{4,t} = CRcells(:,:);
        Percentages(t) = sum(all(cellfun(@(x) x.persistentActivity,RUNS{1,t}))) * 100 / size(RUNS{1,t},2);
        fprintf('Saving results... ');
        %     save('batch_Yiota.mat');
        fprintf('DONE!\n');
    end % for t
    
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
ru = 1 ;% what run of this cluster

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

for ru = 1:10 % what run of this cluster
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
    FFPC(:,:,ru) = TMP./4  ;
end

FFPCtmp = FFPC;

ffstimclu = [];
for i=1:10
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

%% Cross correlation:
close all;
PerClust={};
for ru = 1:10 % what run of this cluster
    for stimclu =1:NC(Sid)
        for cl=1:NC(Sid)
            for c=1:nPC
                if(RUNS{1,stimclu}{c,ru}.clusterID == cl)
                    if(RUNS{1,stimclu}{c,ru}.persistentActivity)
%                         PAPerClust(cl,stimclu) = PAPerClust(cl,stimclu) +1;
                        PerClust{cl,stimclu,ru} = RUNS{1,stimclu}{c,ru}.mv(10001:end);
                    end
                end
            end
        end
    end
end


% gather stimulated cluster membrane depolarisation:
diagonal = {};
for i=1:10
    for d=1:NC(Sid)
        diagonal(d) = PerClust(d,d,i);
    end
    tmp = diagonal( cellfun(@(x) ~isempty(x),diagonal) )  ;
    diagonalALL(i,1:length(tmp)) = tmp  ;
end

diagonalALL = diagonalALL(:)  ;
diagonalALL = diagonalALL(cellfun(@(x) ~isempty(x),diagonalALL))  ;

% Gather all the other culster (non stimulated) depolarisations:
PerClusttmp = PerClust;
for i=1:10
    for d=1:NC(Sid)
        PerClusttmp{d,d,i} = [];
    end
end

PerClusttmp = PerClusttmp(cellfun(@(x) ~isempty(x),PerClusttmp))  ;


% for all the members of diagonalALL make cross correlation with every
% other response.


% Init SPIKY:
para.tmin = 0;
para.tmax = 6000;
para.dts =  10;
para.dtm = 10;
para.select_measures=[0 1     0 0     0 0     0 0];            % Select order of measures

% feed SPIKY
results = SPIKY_no_plot_f_distances_MEX(spikesDend(1,:),para);





% non stimulated
mean(FFPCtmp)
std(FFPCtmp)

% in stimulated cluster:
mean(ffstimclu)
std(ffstimclu)

%%
scatter3(PCsomata(:,1),PCsomata(:,2), PCsomata(:,3), 20, [0,0,1], 'fill', 'o');
axis equal; hold on;