% Built a structured and random network.  
close all; clear all; clc
run.sn = 3;
rng('default') 
rng(run.sn)

% Number of neurons
nPC = 700;
nPV = round(nPC*25/75);% PV mporoun na einai kai 15% !!!
if(nPV==0)
    nPV = 1;
end
nAll = sum([nPC,nPV]);

%  ---- 3D positions ----
inverseDensity =  nthroot(nPC / (75 / 8),3)
cube_dimensions = inverseDensity * 100;
% 3-d points in space. Arguments: #of cells, max seperation distance.
[PCsomata, distPC2PC, distPC2PCwrapped]= CreateRandomNetwork(nPC, cube_dimensions);
distPC2PC = triu(distPC2PC,1) + tril(distPC2PC',-1);
distPC2PCb = distPC2PC;
distPC2PC = triu(distPC2PCwrapped,1) + tril(distPC2PCwrapped',-1);
 
% 3d points for PV. Returns 3d points and distances from PCcells  
% Nassi -Max seperation Distance as PC-PC
[PVsomata, distPV2PC, distPV2PCwrapped] = CreateCubeNetworkPV(cube_dimensions, nPV, PCsomata);
distPV2PCb = distPV2PC;
distPV2PC = distPV2PCwrapped;


run.nPC=nPC;
run.nPV=nPV;
run.nAll=nAll; 
run.stepsperms = 10;



% visualize differences after wrapping (Also do it after clustering?)
tmprn = 0:5:500;
figure;hold on;
title('PC2PC');
plot(tmprn,histc(distPC2PC(logical(triu(ones(nPC),1))),tmprn) / ((nPC^2-nPC)/2),'r');
plot(tmprn,histc(distPC2PCb(logical(triu(ones(nPC),1))),tmprn) / ((nPC^2-nPC)/2),'b');
legend({'No wrapping','Wrapping'});
xlabel('Intersomatic distance (um)');ylabel('Relative frequency');
figure;hold on;
title('PV2PC');
plot(tmprn,histc(distPV2PC(:),tmprn) / (nPC*nPV),'r')
plot(tmprn,histc(distPV2PCb(:),tmprn) / (nPC*nPV),'b')
legend({'No wrapping','Wrapping'});
xlabel('Intersomatic distance (um)');ylabel('Relative frequency');

figure;
scatter3(PCsomata(:,1),PCsomata(:,2),PCsomata(:,3)); 
hold on;
scatter3(PVsomata(:,1),PVsomata(:,2),PVsomata(:,3), 'r');

%% Connect PV-PC neurons 

connProbsPV2PC= @(x) ( 0.4 .* exp(-0.003 * x));
ConnMatPV2PC = connectPV2PC(distPV2PC,connProbsPV2PC); %distance p
ConnMatPC2PV = connectPC2PV(distPV2PC'); % uniform p

% Validate distance dependence of PV 2 PC connections:
rang = 1:5:1000;
rangHisto = histc(distPV2PC(:),rang);
beforeHisto = histc(ConnMatPV2PC(:) .* distPV2PC(:),rang) ./ rangHisto;

% Test if the PC2PV reciprocals are 20% of all pairs. 
% Based on Otsuka, Kaywagushi-Figure2B. Rearrange connections accordingly.
ctr = 80000;
while ctr && (0.1 < abs(((sum(sum((ConnMatPC2PV & ConnMatPV2PC')))*100) / numel(ConnMatPC2PV)) - 20))
    ctr
    % if reciprocals are more/less, move connections to reach 20%:
    foundIt = 0;
    while ~foundIt
        ip = ceil(rand(1)*size(ConnMatPC2PV,1));
        jp = ceil(rand(1)*size(ConnMatPC2PV,2));
%         If there PC-to-PV connection but not PV to PC
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

% determine no of PV to PC connections needed to meet the observations:
diffFronObserved = beforeHisto' - connProbsPV2PC(rang);
candidates = find(diffFronObserved < 0) ;
[S, I] = sort(diffFronObserved(candidates));
connsNeeded = abs(ceil(diffFronObserved(candidates(I)) .* rangHisto(candidates(I))' ));
for k = 1:length(candidates)
    rangeBin = rang(candidates(I(k)):candidates(I(k))+1);
    idx = find( (rangeBin(1) <= distPV2PC(:)) & (distPV2PC(:) < rangeBin(2)) );
    while (connsNeeded(k) > 0) && (sum(~ConnMatPV2PC(idx)))
            tmp = find(~ConnMatPV2PC(idx));
            ConnMatPV2PC(idx(tmp(randperm(length(tmp),1)))) = 1;
            connsNeeded(k) = connsNeeded(k) - 1 ;
    end
end

figure;hold on;
plot(beforeHisto);
plot(connProbsPV2PC(rang),'r');
afterHisto = histc(ConnMatPV2PC(:) .* distPV2PC(:),rang) ./ rangHisto;
plot(afterHisto,'k');

% Otsuka: No connection 0.7, reciprocal 0.2, PC-IN 3%, IN-PC 7%
Con_prob(1,1)=(sum(sum((ConnMatPC2PV==0 & ConnMatPV2PC'==0)))*100) / numel(ConnMatPC2PV);
Con_prob(2,1)=(sum(sum((ConnMatPC2PV==1 & ConnMatPV2PC'==1)))*100) / numel(ConnMatPC2PV);
Con_prob(3,1)=(sum(sum((ConnMatPC2PV==1 & ConnMatPV2PC'==0)))*100) / numel(ConnMatPC2PV);
Con_prob(4,1)=(sum(sum((ConnMatPC2PV==0 & ConnMatPV2PC'==1)))*100) / numel(ConnMatPC2PV);

cm = [1,1,1;0,0,0;0.2,0.2,0.2;0.5,0.5,0.5];
figure;bar([Con_prob,Con_prob]','stacked');
xlabel('PC to PV');ylabel('Percentage');
colormap(cm);

%PV-PV. Connection Probability: Galarreta, Hestrin, 1999 
% (somatosensory/visual L5):
% 66% electrical coupling <=80? appart.
% 18% GABAergic connections
% 4.5% Reciprocal GABAergic connections
[ConnMatPV2PV, gapmat]=connectPV2PV(nPV);
%Label the gap junctions for NEURON:
start = 2;
ctr = 1;
for k = 1:nPV
    for j = start:nPV
        if gapmat(k,j)
            gapmat(k,j) = ctr;
            ctr = ctr + 2;
        end
    end
end
gtmp = gapmat';
gapmat(logical(tril(ones(nPV),-1))) = gtmp(logical(tril(ones(nPV),-1)));

% Set indices for ease of mind:
pc = size(PCsomata,1);
pv = size(PCsomata,1) + size(PVsomata,1);
   
% Populate the final all-to-all connectivity matrix. In this matrix rows
% are the sources, columns are the targets.
AllConnMat = zeros(nAll);
% PCs to PVs
AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
% PVs to PVs + autapses 
AllConnMat(pc+1:pv,pc+1:pv) = connectPV2PV(nPV);
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;

%%
total_prob = @(x) 0.22.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x))); % REMOVE THE MULTIPLICATION FACTOR!!
% SOS remove the multiplication factor because of error in the paper: SOS
upto = 400;
fact = 1;%665;
figure;hold on;
plot(total_prob(0:upto)*fact,'r');
plot(connProbsPC2PC(0:upto)*fact,'k');
plot(recipProbsPC2PC(0:upto)*fact,'b');
legend({'Overall','Unidirectional','Reciprocal'});
xlabel('Intersomatic distance (um)');ylabel('Connection probability');

%% Create Structured Network:
myrange = 1:10:max(distPC2PC(:))+10;
orderedDist = triu(distPC2PC,1);
orderedDist = orderedDist(logical(triu(ones(nPC,nPC),1)));
distBins = histc( orderedDist, myrange);
% Create Pd:
Pd = zeros(nPC);
for i=1:nPC
    for j=1:nPC
        if (i~=j)
            Pd(i,j)=total_prob(distPC2PC(i,j));
        else
            Pd(i,j)=0;
        end
    end
end

% Verify complex connectivity:
E = zeros(nPC);
start = 2;
for i=1:nPC
    for j=start:nPC
        xx=distPC2PC(i,j);
        pnr = connProbsPC2PC(xx);
        pr  = recipProbsPC2PC(xx);
        randomNum = rand(1);
        if(randomNum < pr)
            E(i,j) = 1;
            E(j,i) = 1;
        elseif ((randomNum >= pr) && (randomNum < pr+(pnr/2)))
            E(i,j) = 1;
        elseif ((randomNum >= pr+(pnr/2)) && (randomNum < (pr+pnr)))
            E(j,i) = 1;
        end
    end
    start = start+1;
end
Er = E|E';
Er = Er & logical(triu(ones(size(Er)),1));
sum(Er(:)) / ((nPC^2-nPC)/2)
E = logical(E);

% run iterative rearrangement: 
iters = 300;
Es = false(nPC,nPC,iters);
F_prob = zeros(nPC,nPC);
CN_prob = zeros(nPC,nPC);
Es(:,:,1) = E;


sumPd = sum(sum(Pd));
for t=2:iters
    t
    NCN= m_commonNeighbors(Es(:,:,t-1));
    CN_prob=((NCN)./max(NCN(~eye(nPC))));
    
    F_prob = CN_prob .* Pd; 
    F_prob = F_prob * ( sumPd / sum(F_prob(:)) );
    
    xx=(log(0.22)-log(F_prob))/0.0052;
    xx(xx<0) = 0;
    pr= recipProbsPC2PC(xx);
    pnr=connProbsPC2PC(xx);
    randomNum = rand(nPC);
    tmpRecip = (randomNum<pr) & triu(ones(nPC),1);
    tmpRecip  = tmpRecip | tmpRecip';
    tmpConn = ((pr <= randomNum) & (randomNum < pr+(pnr)));
    tmpConn = tmpConn & triu(ones(nPC),1);
    
    tmpV = find(tmpConn);
    tmpI = false(1,length(tmpV));
    tmpI(rand(1,length(tmpV))<0.5) = 1;
    tmpConnU = false(nPC);
    tmpConnL = false(nPC);
    tmpConnU(tmpV(tmpI)) = 1;
    tmpConnL(tmpV(~tmpI)) = 1;
    tmpConnL = tmpConnL';
    Es(:,:,t) = tmpConnU | tmpConnL | tmpRecip;
end

nUniqueUnorderedPairs = (nPC^2 - nPC)/2;
for t=1:iters
    t
    connected = Es(:,:,t)|Es(:,:,t)';
    reciprocal = Es(:,:,t)&Es(:,:,t)';
    connected = connected & logical(triu(ones(size(connected)),1));
    
    connectedUnorderedPairsPercentage(t) = sum(connected(:)) / nUniqueUnorderedPairs;% overall connection probability
    reciprocalUnorderedPairsPercentage(t) = (sum(reciprocal(:))/2) / nUniqueUnorderedPairs;
    
    [~,cc(t),~] = clust_coeff(Es(:,:,t));

end

figure;plot(connectedUnorderedPairsPercentage);
ylabel('Overall connection probability');xlabel('Iterations');
figure;plot(reciprocalUnorderedPairsPercentage);
ylabel('Reciprocal connection probability');xlabel('Iterations');
figure;plot(cc);
ylabel('Average local clusterin coefficient');xlabel('Iterations');


PC2PC_str = Es(:,:,end);
prob_conn_rnd_ind = sum(PC2PC_str(:)) / (nPC^2-nPC);
[Ernd] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind');
PC2PC_d = E;
PC2PC_rnd = Ernd;

% Average connection probability:
sum(sum( (PC2PC_rnd|PC2PC_rnd') & logical(triu(ones(size(PC2PC_rnd)),1)) )) / ((nPC^2-nPC)/2)
sum(sum( (PC2PC_str|PC2PC_str') & logical(triu(ones(size(PC2PC_str)),1)) )) / ((nPC^2-nPC)/2)

% Get total configurations:
run.configuration_rnd = AllConnMat;
run.configuration_rnd(1:run.nPC, 1:run.nPC) = PC2PC_rnd;
for k = 1: run.nPC
    run.configuration_rnd(k,k) = 1;
end
run.configuration_str = AllConnMat;
run.configuration_str(1:run.nPC, 1:run.nPC) = PC2PC_str;
for k = 1: run.nPC
    run.configuration_str(k,k) = 1;
end
%% Create weight matrices:

% TODO:
% Generate weights matrices:
% Kayaguchi : more potent synapses in reciprocal pairs
% Perin: more potent synapses with more clustering coeff:
% We try to combine them

% Extrapolating Perin's findings, more CN, more shift in distribution:
snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
rang = 0:0.01:8;
maxNeighbors = max(max([m_commonNeighbors(PC2PC_d),m_commonNeighbors(PC2PC_rnd),m_commonNeighbors(PC2PC_str)]));
alpha = 3*ones(1,maxNeighbors);%linspace(3,0.1,maxNeighbors)%.*exp(linspace(0,0.5,maxNeighbors)).^-2;%[20,10,2,1];
maxSigma = 5 ;
minSigma = 0.1 ;
expClimb = (exp(linspace(0,0.5,maxNeighbors)).^5)-1;
expClimb = (expClimb/max(expClimb)*maxSigma)+minSigma;
m = ones(1,maxNeighbors).*expClimb.*2;%exp(linspace(0,0.6,maxNeighbors))-1;
m = m-m(1);
sigma = 0.5*ones(1,maxNeighbors).*expClimb%.*(1./(1+exp(-linspace(-2,5,maxNeighbors))))%.*exp(linspace(0,0.5,maxNeighbors)).^-10;%[0.15, 0.2, 4, 6];
normFactors = ones(1,maxNeighbors);%[0.0045, 0.006, 0.02, 0.02];
% LogNormal PDF parameters:
parmhat(1,1:maxNeighbors) = linspace(0.1,10,maxNeighbors);
parmhat(2,1:maxNeighbors) = 2.6;
CNsnormal = [];
% Probabilities changed to replicate mean EPSP amplitude per connections in
% cluster (Perin et al, Fig.6A):
maxCDF = zeros(maxNeighbors,length(rang));
for k=1:maxNeighbors
    CNsnormal(k,:) = snormal(rang, m(k), alpha(k), sigma(k))*normFactors(k);
    maxCDF(k,:) = cumsum(CNsnormal(k,:))/max(cumsum(CNsnormal(k,:)));
end
figure;
plot(CNsnormal(:,:)');hold on;
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));

CN_str = m_commonNeighbors(PC2PC_str);
weights_str = ones(nPC,nPC);
tmp_weights_str = ones(nPC,nPC);
for k=1:maxNeighbors
    k
    % Even more efficient: get them instantly from maxCDF:
    randomSampling = rand(1,sum(sum(CN_str==k-1)));
    tmp = histc(randomSampling,[0,maxCDF(k,:)]);
    idxvec = zeros(1,sum(tmp));
    cumvec = [0,cumsum(tmp)];
    for kk=1:length(cumvec)-1
        idxvec(cumvec(kk)+1:cumvec(kk+1)) = kk;
    end
    tmp_weights_str(CN_str==k-1) = rang(idxvec(randperm(length(idxvec))));
    
end
weights_str(:,:) = tmp_weights_str .* PC2PC_str;

% cluster and then poso einai to mean weight per cluster?
nweights_str = (weights_str/max(weights_str(:)))*5;
run.weights_str = nweights_str;
% Random configuration weights:
% Structured configuration weights distribution:
tmp = nweights_str .* PC2PC_str;
tmp(tmp==0) = [];
% histcounts(tmp)
run.weights_rnd = PC2PC_rnd .* ones(run.nPC) * median(tmp);
%% Run Clustering:

mr = -200:5:-60; % preference range
tmpidx = cell(1,length(mr));
dpsim = cell(1,length(mr));
[X,Y] = meshgrid(1:nPC);
similarity_str = zeros((nPC^2-nPC)/2,3);
similarity_str(:,1) = X(logical(triu(ones(nPC),1)));
similarity_str(:,2) = Y(logical(triu(ones(nPC),1)));
CN_str = m_commonNeighbors(PC2PC_str);
similarity_str(:,3) = CN_str(logical(triu(ones(nPC),1)));
for k=1:length(mr)
    mr(k)
[tmpidx{k},~,dpsim{k},~]=apclusterSparse(similarity_str,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
% [targ_str,~,labels_str] = unique(tmpidx{k}(:,end));

% Random network:
CN_rnd = m_commonNeighbors(PC2PC_rnd);
tmpidx_rnd = cell(1,length(mr));
dpsim_rnd = cell(1,length(mr));
similarity_rnd = zeros((nPC^2-nPC)/2,3);
similarity_rnd(:,1) = X(logical(triu(ones(nPC),1)));
similarity_rnd(:,2) = Y(logical(triu(ones(nPC),1)));
similarity_rnd(:,3) = CN_rnd(logical(triu(ones(nPC),1)));
for k=1:length(mr)
    mr(k)
[tmpidx_rnd{k},~,dpsim_rnd{k},~]=apclusterSparse(similarity_rnd,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
% [targ_rnd,~,labels_rnd] = unique(tmpidx_rnd{k}(:,end));

% figure;plot(cellfun(@(x) x(end), dpsim_rnd))
figure;hold on;
plot(cellfun(@(x) length(unique(x(:,end))),tmpidx_rnd ))
plot(cellfun(@(x) length(unique(x(:,end))),tmpidx ))
legend({'RND','STR'})

pref_str = find(cellfun(@(x) length(unique(x(:,end))),tmpidx )==7);
pref_str = pref_str(1);

pref_rnd = find(cellfun(@(x) length(unique(x(:,end))),tmpidx_rnd )==7);
pref_rnd = pref_rnd(1);

[~,~,clusterLabels_rnd]=unique(tmpidx_rnd{pref_rnd}(:,end))
[~,~,clusterLabels_str]=unique(tmpidx{pref_str}(:,end))

for k=1:7
   nc_rnd(k) =  sum(clusterLabels_rnd==k);
   nc_str(k) =  sum(clusterLabels_str==k);
end


run.nClusters_rnd = length(nc_rnd);
run.nClusters_str = length(nc_str);
run.nCellsInCluster_rnd = nc_rnd;
run.nCellsInCluster_str = nc_str;

run.stimulatedCells_rnd = cell(1,run.nClusters_rnd);
run.stimulatedCells_str = cell(1,run.nClusters_str);
for k = 1:run.nClusters_rnd
    run.stimulatedCells_rnd{k} = find(clusterLabels_rnd==k);
end
for k = 1:run.nClusters_str
    run.stimulatedCells_str{k} = find(clusterLabels_str==k);
end

%% Save configuration:
% With different serial number
save(sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\NetworkCreation_SN%d.mat',run.sn),'-v7.3');

%% Export to NEURON environment:

pathto = '\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\';
exportNetworkParameters(run, 'str', pathto);
exportNetworkParameters(run, 'rnd', pathto);
exportStimulationParameters(run, pathto);

