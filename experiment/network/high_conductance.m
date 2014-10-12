
aa=9
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

    
    exprun=10;
    stimend = 1500;
    nPC = 75;%225; %216;
%     clstr = 1; % Cluster type
    
% % Create tissue sanbox of given dimensions to fit all the pertitnent
% % neuronal types in it, in correct densities.
% tSandbox.side = 250;
% PVsomata = CreateCubeNetworkPV(250, 12.5);
%     scatter3(PVsomata(:,1),PVsomata(:,2),PVsomata(:,3));

    factorOverallProb = 0.7;

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
    
    % Pyramidals connect to all
    PC2PC = zeros(nPC,nPC,2);
    connProbsPC2PC = connProbsPC2PC*factorOverallProb;
    recipProbsPC2PC = recipProbsPC2PC *factorOverallProb;
    NincomingProbs = NincomingProbs *factorOverallProb;
    NoutgoingProbs = NoutgoingProbs *factorOverallProb;
    PC2PC(:,:,1) = create_graph_CN(distPC2PC,connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    
    
    prob_conn_rnd_joint = sum(sum(tril(PC2PC(:,:,1) | PC2PC(:,:,1)') )) / ((numel(PC2PC(:,:,1))-nPC)/2); % individual probability
    prob_conn_rnd_ind = sum(sum(PC2PC(:,:,1))) / (numel(PC2PC(:,:,1))-nPC); % joint probability
    
    PC2PC(:,:,2) = create_graph_WS(distPC2PC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
    
    [mcccs,~,cccs] = clust_coeff( PC2PC(:,:,1) );
    [mcccr,~,cccr] = clust_coeff( PC2PC(:,:,2) );
    [mcccr,mcccs]
    
    
    MAXPC2PC = [];
    CCDIFF = [];
    CCMAX = 0;
    
    for i=1:100
        i
        PC2PC(:,:,1) = create_graph_CN(distPC2PC,connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
        prob_conn_rnd_ind = sum(sum(PC2PC(:,:,1))) / (numel(PC2PC(:,:,1))-nPC); % joint probability
        PC2PC(:,:,2) = create_graph_WS(distPC2PC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
        [mcccs,~,cccs] = clust_coeff( PC2PC(:,:,1) );
        [mcccr,~,cccr] = clust_coeff( PC2PC(:,:,2) );
        CCDIFF(i) = mcccs - mcccr
        if(CCMAX < CCDIFF(i))
            disp(sprintf('Replacing %f with %f',CCMAX,CCDIFF(i)))
            CCMAX = CCDIFF(i)
            MAXPC2PC(:,:,1) = PC2PC(:,:,1); 
            MAXPC2PC(:,:,2) = PC2PC(:,:,2); 
        end
    end
    max(CCDIFF)
    PC2PC(:,:,1) = MAXPC2PC(:,:,1); 
    PC2PC(:,:,2) = MAXPC2PC(:,:,2); 
    
    
    sum(sum(PC2PC(:,:,1)))
    sum(sum(PC2PC(:,:,2)))
    
    
    sum(sum(PC2PC(:,:,1) & PC2PC(:,:,1)'))
    sum(sum(PC2PC(:,:,2) & PC2PC(:,:,2)'))

    
    % ---- affinity propagation similarity measure: ----
    
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
    AllConnMat(1:pc,pv+1:end) = rand(nPC, nAllCells-pv) > 0.5; % Connected to all interneurons: too much. maybe less connectivity?
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
    
    
    % Generate stimulation data:
    nIncoming=[200,200];
    randShift = 0.1;
    nStim = 30; % how many Hz of stimulation?
    tstop = 1000;
    spikesDend = {};
    spikesApic = {};
    for c=1:pc
        i=1;
        while i <= nIncoming(1)
            spikesDend(c,i)={linspace(0,tstop,nStim)};
            spikesDend{c,i} = round(spikesDend{c,i} + ((rand(1,nStim)-0.5)*randShift));
            spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}>0);
            spikesDend{c,i} = spikesDend{c,i}(spikesDend{c,i}<tstop);
            spikesDend{c,i} = sort((spikesDend{c,i}),'ascend');
            if(length(unique(spikesDend{c,i})) ~= length(spikesDend{c,i}) )
                continue;
            end
            i = i +1;
        end
        i=1;
        while i <= nIncoming(2)
            spikesApic(c,i)={linspace(0,tstop,nStim)};
            spikesApic{c,i} = round(spikesApic{c,i} + ((rand(1,nStim)-0.5)*randShift));
            spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}>0);
            spikesApic{c,i} = spikesApic{c,i}(spikesApic{c,i}<tstop);
            spikesApic{c,i} = sort((spikesApic{c,i}),'ascend');
            if(length(unique(spikesApic{c,i})) ~= length(spikesApic{c,i}) )
                continue;
            end
            i = i +1;
        end
    end
    
% NASSI 9Neurons SYNCHRONOUS




 % Background activity must result in ~1.2Hz in the post-synaptic cell as in
 % Carl C.H.Petersen Neuron, 2013
tstop = 5000;
many_times = 1;
num_exc = 140;%100;
num_inha = 75;%75
num_inhb = 30;%30
num_inhc = 30;%30

spikes_basal={};
spikes_proximal={};
spikes_apical={};
spikes_interneurona={};
spikes_interneuronb={};
spikes_interneuronc={};

% for many_times=1:runs+1
%     rng(many_times)
    stimPoisson = find(poissrnd(0.001,tstop,1))  ;%0.006 mporei na einai poly gia to structured
%     stimPoisson = sort([stimPoisson,rand(tstop]);
% stimPoisson = [200:400:3000]' ;
TID=100;

% POUSTIA terastiwn diastasewn:
% Apo to Poisson trace bgazw ta events pou einai pio konta apo time TID
stimPoisson(diff(stimPoisson)<(TID*2)) = [];

% stimPoisson = [TID:TID:tstop]' ;
    
% Mipws oi synapseis pou dinw gia background na einai anti-clustered?
    for c=1:nPC
       poil = length(stimPoisson) ;
        for i = 1:num_exc
            tmp = sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ;
            spikes_basal(many_times,c,i) = {unique(tmp(tmp<tstop))};
            spikes_basal{many_times,c,i}(spikes_basal{many_times,c,i}<0) = [];
        end
        for i = 1:num_exc
            tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
            spikes_proximal(many_times,c,i) = {unique(tmp(tmp<tstop))};
            spikes_proximal{many_times,c,i}(spikes_proximal{many_times,c,i}<0) = [];
        end
        for i = 1:num_exc
            tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
            spikes_apical(many_times,c,i) = {unique(tmp(tmp<tstop))};
            spikes_apical{many_times,c,i}(spikes_apical{many_times,c,i}<0) = [];
        end
    end
    
    for c = 1:nPV
        for i = 1:num_inha
            tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
            spikes_interneurona(many_times,c,i) =  {unique(tmp(tmp<tstop))};
            spikes_interneurona{many_times,c,i}(spikes_interneurona{many_times,c,i}<0) = [];
        end
    end
    
    for c = 1:nCB
        for i = 1:num_inhb
            tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
            spikes_interneuronb(many_times,c,i) =  {unique(tmp(tmp<tstop))};
            spikes_interneuronb{many_times,c,i}(spikes_interneuronb{many_times,c,i}<0) = [];
        end
    end
    
    for c = 1:nCR
        for i = 1:num_inhc
            tmp =  sort([round(((rand(poil,1)-0.0)*TID)) + stimPoisson;rand(tstop/10000,1)*tstop],'ascend') ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
            spikes_interneuronc(many_times,c,i) =  {unique(tmp(tmp<tstop))};
            spikes_interneuronc{many_times,c,i}(spikes_interneuronc{many_times,c,i}<0) = [];
        end
    end
% end

c=1;figure();
subplot(4,1,1)
for syn = 1:num_exc
        for l=1:size(spikes_basal{many_times,c,syn},1)
            plot([spikes_basal{many_times,c,syn}(l),spikes_basal{many_times,c,syn}(l)],[syn-1,syn],'k');hold on;
        end  
end
ylabel('Excitatory')
set(gca,'ytick',[])

subplot(4,1,2)
for syn = 1:num_inha
        for l=1:size(spikes_interneurona{many_times,c,syn},1)
            plot([spikes_interneurona{many_times,c,syn}(l),spikes_interneurona{many_times,c,syn}(l)],[syn-1,syn],'k');hold on;
        end  
end
ylabel('PV')
set(gca,'ytick',[])

subplot(4,1,3)
for syn = 1:num_inhb
        for l=1:size(spikes_interneuronb{many_times,c,syn},1)
            plot([spikes_interneuronb{many_times,c,syn}(l),spikes_interneuronb{many_times,c,syn}(l)],[syn-1,syn],'k');hold on;
        end  
end
ylabel('CB')
set(gca,'ytick',[])

subplot(4,1,4)
for syn = 1:num_inhc
        for l=1:size(spikes_interneuronb{many_times,c,syn},1)
            plot([spikes_interneuronb{many_times,c,syn}(l),spikes_interneuronb{many_times,c,syn}(l)],[syn-1,syn],'k');hold on;
        end  
end
ylabel('CR')
set(gca,'ytick',[])

figure();
voidvec = [];
for syn = 1:num_exc
    voidvec = [voidvec; spikes_basal{many_times,c,syn}] ;
end
plot( histc(voidvec,1:10:tstop) )

% for syn = 1:i
%     plot([1,2],[3,4])
%     scatter(spikes_basal{many_times,c,syn},ones(length(spikes_basal{many_times,c,syn}),1)*syn,'+');hold on;
% end

exportDetailedBackgroundStimParams(spikes_basal,spikes_proximal, spikes_apical,...
    spikes_interneurona, spikes_interneuronb, spikes_interneuronc, mypath);

%     % GABAa inhibition in pyramidals' somata:
%     % approximate No of PV in range of each pyramidal connected
%     % propabilisticaly with distance as in R.Yuste, 2011
%     PVdummies = rand(nPV*5*5*3,3);  % 226.6 per mm squared (?!?) paper?
%     PVdummies = PVdummies .* repmat([750,750,450],length(PVdummies),1);
%     distPVd2PC = distancePV2PC(PVdummies,[375,375,225]);
%     distPVd2PC(distPVd2PC >= 400) = 0;
%     ConnMatPVd2PC = connectPV2PC(distPVd2PC,connBinsPV2PC,connProbsPV2PC);
%     sum(ConnMatPVd2PC) % No of inhibitory GABAa synapses to a pyramidal

%     % Generate background inhibition data:
%     nIncoming=[396];% hardcoded to avoid run to run variability
%     randShift = 1000/6;
%     nStim = 30; % ~6Hz as in Alain Destexhe & Denis Pare, 1999
%     tstop = 5000;
%     spikesGABAa = {};
%     for c=1:pc
%         i=1;
%         while i <= nIncoming(1)
%             spikesGABAa(c,i)={linspace(0,tstop,nStim)};
%             spikesGABAa{c,i} = round(spikesGABAa{c,i} + ((rand(1,nStim)-0.5)*randShift));
%             spikesGABAa{c,i} = spikesGABAa{c,i}(spikesGABAa{c,i}>0);
%             spikesGABAa{c,i} = spikesGABAa{c,i}(spikesGABAa{c,i}<tstop);
%             spikesGABAa{c,i} = sort((spikesGABAa{c,i}),'ascend');
%             if(length(unique(spikesGABAa{c,i})) ~= length(spikesGABAa{c,i}) )
%                 continue;
%             end
%             i = i +1;
%         end
%     end
    
    % run for the structured network and for the random network as
    % well:
    
    RUNS_str = {};
    RUNS_rnd = {};
    
%     export network parameters for that experiment:
    cd(mypath);
    system(sprintf('mkdir experiment_%d',aa));
    exportNetworkPositions(PCsomata,PVsomata,sprintf('%sexperiment_%d/',mypath,aa));
    % Export pyramidal's clusters ID:
%     exportNetworkCluster([idx_str,idx_rnd],sprintf('%sexperiment_%d/',mypath,aa));
    % Export stimulation for each pyramidal:
    exportNetworkStimulation(spikesDend,spikesApic,sprintf('%sexperiment_%d/',mypath,aa));
%     exportNetworkInhibition(spikesGABAa,sprintf('%sexperiment_%d/',mypath,aa));
    
%     save(sprintf('%sexperiment_%d/experiment.mat',mypath,aa),'-v7.3');
    
    
%     cellsPerCluster_str = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
%     cellsPerCluster_rnd = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
%     
%     stimCellsPerCluster = min(min(cellsPerCluster_str{:}),min(cellsPerCluster_rnd{:}));
    
StimVect_str = zeros(nPC,1);
% stimulation parameters for the test runs:
% PcellStimList.x[14] = 1
% PcellStimList.x[20] = 1
% PcellStimList.x[30] = 1
% PcellStimList.x[32] = 1
% PcellStimList.x[37] = 1
% PcellStimList.x[57] = 1
% PcellStimList.x[66] = 1
% StimVect_str([15,21,31,33,38,58,67],1) = 1;
% Random stimulation that gives no PA, but PA in STR:
% StimVect_rndA = StimVect_str;
% AllConnMat_strA = AllConnMat_str;
% AllConnMat_rndA = AllConnMat_rnd;
StimVect_str([round(rand(5,1)*75)],1) = 1;

AllConnMat_str = AllConnMat;
AllConnMat_str(1:pc,1:pc) = PC2PC(:,:,1);
AllConnMat_str(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1; 

AllConnMat_rnd = AllConnMat;
AllConnMat_rnd(1:pc,1:pc) = PC2PC(:,:,2);
AllConnMat_rnd(sub2ind([nPC,nPC],[1:pc],[1:pc])) = 1;

AllWeightMat_str = AllWeightMat;


cellsPerCluster_str = {histc(labels_str(:,Sid),1:NC_str(Sid))'};
cellsPerCluster_rnd = {histc(labels_rnd(:,Sid),1:NC_rnd(Sid))'};
stimCellsPerCluster = min(min(cellsPerCluster_str{:}),min(cellsPerCluster_rnd{:}));
% Isolate one cluster:
t=1
[r,~,~] = find(labels_str(:,Sid) == t);
StimVect_str = zeros(nPC,1);
tmpr = randperm(length(r));
StimVect_str(r(tmpr(1:stimCellsPerCluster))) = 1;  % 6==NUmber of total cells stimulated in each cluster.


[r,~,~] = find(labels_rnd(:,Sid) == t);
StimVect_rnd = zeros(nPC,1);
tmpr = randperm(length(r));
StimVect_rnd(r(tmpr(1:stimCellsPerCluster))) = 1;
    


%         AllWeightMat_str = AllWeightMat_str .* AllConnMat_str ;
AllWeightMat_str = ones(size(AllConnMat_str));

%  ---- Export parameters to NEURON ----
AllConnMat_str_TEST = AllConnMat_str;
AllConnMat_str_TEST(pc+1:pv,1:pc) = ones(size(ConnMatPV2PC));
%         AllConnMat_str_TEST(r(tmpr(1:6)),r(tmpr(1:6))) = 1; 
%         AllConnMat_str_TEST(1:end,35) = 0; 
exportNetworkParameters(AllCells,AllConnMat_rnd,AllWeightMat_str,sprintf('%sexperiment_%d/',mypath,aa));
exportStimulationParameters(AllCells,StimVect_rnd,sprintf('%sexperiment_%d/',mypath,aa));

t=1; % Initial without vclamps
VCLAMP = 0;

fprintf('Running NEURON... ');
unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
tic;
%         [nrnStatus,nrnCmdOut{t}] = unix('./parallel_pfc');
EXPERIMENT = 1;  % 0 = random, 1= structured gia to $#$#% NEURON
%         [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpiexec -np 8 %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" final.hoc',mypath2,t,EXPERIMENT,aa));
[nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" -c "VCLAMP=%f" final.hoc',mypath2,t,EXPERIMENT,aa, VCLAMP));
runtime = toc

if( ~findstr('Success',nrnCmdOut{t}))
    continue;
end

RESTING = -66;
test_vclamp = -10:-10:-70;

for t=2:2
    fprintf('Running NEURON... ');
    unix(sprintf('kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk ''{print$2}''`'));
    tic;
    EXPERIMENT = 1;  % 0 = random, 1= structured gia to $#$#% NEURON
    VCLAMP = RESTING + test_vclamp(t-1);
    [nrnStatus,nrnCmdOut{t}] = unix(sprintf('mpirun -n 24 -mca btl ^openib %smechanism_complex/x86_64/special -mpi -nobanner -c "PARALLEL=1" -c "SIMPLIFIED=0" -c "CLUSTER_ID=%d" -c "EXPERIMENT=%d"  -c "AA=%d" -c "VCLAMP=%f" final.hoc',mypath2,t,EXPERIMENT,aa, VCLAMP));
    runtime = toc
    
    if( ~findstr('Success',nrnCmdOut{t}))
        continue;
    end
    
end

fprintf('DONE!\n');

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
%choose structured or random:
% expString = 'RND';
expString = 'STR';

PCcells_str = {};
PCcells_str_spk = {};
PCcells_str_i = {};
PCcells_str_dist= {};
for t=1:1
    for ru = 1:exprun
        for c=1:pc
            mycell = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d.txt',mypath,aa,expString,t,c-1,ru-1)),10);
            mycell.position = PCsomata(c,1:3);
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
        if t>1
            PCcells_str_i{t-1,ru} = ncell(load(sprintf('%sexperiment_%d/%s/%d_%d_%d_i.txt',mypath,aa,expString,t,0,ru-1)),10);
        end
    end
    RUNS_str{1,t} = PCcells_str(:,:);
end


fprintf('DONE!\n');

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

