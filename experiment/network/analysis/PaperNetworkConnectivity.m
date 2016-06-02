% Stefanos, 2014 
% Built a structured and random network.  
close all; clear all; clc
% Random seed
aa=1;
rng('default') 
rng(aa)

%Set paths
% mypath = '/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files';
% mypath2 = '/home/cluster/papoutsi/Desktop/STEFANOS/nassi/';
% addpath('/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files/adapt_apV3')
% addpath('/home/cluster/papoutsi/Desktop/STEFANOS/nassi/experiment/network/matlab_files/code')

% Number of neurons
nPC = 700;
nPV = round(nPC*25/75);% PV mporoun na einai kai 15% !!!
if(nPV==0)
    nPV = 1;
end
nAllCells = sum([nPC,nPV]);
AllCells = [nPC , nPV , nAllCells];

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
[PVsomata, distPV2PC, distPC2PCwrapped] = CreateCubeNetworkPV(cube_dimensions, nPV, PCsomata);
distPV2PCb = distPV2PC;
distPV2PC = distPC2PCwrapped;

% visualize differences after wrapping (Also do it after clustering?)
tmprn = 0:5:500;
figure;hold on;
title('PC2PC');
plot(tmprn,histc(distPC2PCb(logical(triu(ones(N),1))),tmprn) / ((N^2-N)/2),'b');
plot(tmprn,histc(distPC2PCwrapped(logical(triu(ones(N),1))),tmprn) / ((N^2-N)/2),'r');
figure;hold on;
title('PV2PC');
plot(tmprn,histc(distPV2PC(:),tmprn) / (nPC*nPV),'b')
plot(tmprn,histc(distPV2PCwrapped(:),tmprn) / (nPC*nPV),'r')

figure;
scatter3(PCsomata(:,1),PCsomata(:,2),PCsomata(:,3)); 
hold on;
scatter3(PVsomata(:,1),PVsomata(:,2),PVsomata(:,3), 'r');

% optional: for big enough networks wrap the dimensions and check what
% changes. 



%% Connect PV-PC neurons 
% Distance of each interneuron from Pyramidal and connect based on
% probability from Yuste 2011:    
% gia to connectivity twn PV 2 PC o Packer et al., 2011:
% elegxei gia ena region gyrw apo to ka8e PC. Epomenws ta histograms einai
% swsta gia ta PV2PC (ka8ws briskontai sto idio layer). Oso gia ta
% CB(SOM)2PC pou briskontai se diaforetika layers mporoume na eikasoume
% oti:
% * Apo to sxhma twn SOM (Packer et al., 2013), oso metakineisai sto 
% transverse layer oi tuft dendrites enws PC 8a exoun tis idies pi8anotites
% gia overlap opos exoun kai ta PC pou briskontai sto idio layer. Epomenws
% metraw san distance CB2PC mono to distance sto transverse plane kai
% xrisimopoiw tis pi8anotites pou dinei o Packer et al., 2013, Fig4D
% Packer, A. M., & Yuste, R. (2011). Dense, unspecific connectivity of 
% neocortical parvalbumin-positive interneurons: a canonical microcircuit 
% for inhibition? The Journal of Neuroscience : The Official Journal of the
% Society for Neuroscience, 31(37), 1326071. 
% doi:10.1523/JNEUROSCI.3131-11.2011

% Packer, A., McConnell, D., Fino, E., & Yuste, R. (2013). Axo-dendritic 
% overlap and laminar projection can explain interneuron connectivity to 
% pyramidal cells. Cerebral Cortex, (December), 27902802. 
% doi:10.1093/cercor/bhs210

% Connection of PCs to PVs as in :
% Cortical inhibitory cell types differentially form intralaminar
% and interlaminar subnetworks with excitatory neurons, Otsuka Takeshi
% Kawaguchi Yasuo, 2009
% Random, with probability 23%
  
% PV-PC
% Nassi
connProbsPV2PC= @(x) ( 0.4 .* exp(-0.003 * x));
% na ypologisw kai to distance PV2PC me wrapping and see what happens
ConnMatPV2PC = connectPV2PC(distPV2PC,connProbsPV2PC); %distance p
ConnMatPC2PV = connectPC2PV(distPV2PC'); % uniform p

% Validate distance dependence of PV 2 PC connections:
rang = 1:5:1000;
rangHisto = histc(distPV2PC(:),rang);
beforeHisto = histc(ConnMatPV2PC(:) .* distPV2PC(:),rang) ./ rangHisto;



% Test if the PC2PV reciprocals are 20% of all pairs. Based on Otsuka, Kaywagushi-Figure2B. Rearrange connection accordingly.:
ctr = 80000;
while ctr && (0.1 < abs(((sum(sum((ConnMatPC2PV & ConnMatPV2PC')))*100) / numel(ConnMatPC2PV)) - 20))
    
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

diffFronObserved = beforeHisto' - connProbsPV2PC(rang);
candidates = find(diffFronObserved < 0) ;
% tmp = diffFronObserved;
% tmp(~candidates) = 0;
[S, I] = sort(diffFronObserved(candidates));
% determine no of connections needed to meet the observations:
% connsNedded = ((numel(ConnMatPC2PV) * 7 ) / 100) - sum(sum((ConnMatPC2PV==0 & ConnMatPV2PC'==1))) ;
% tmp = histc(ConnMatPV2PC(:) .* distPV2PC(:),rang)';
% connsNedded = ceil(connProbsPV2PC(candidates(I)) .* rangHisto(candidates(I))' ) - tmp(I);

connsNedded = abs(ceil(diffFronObserved(candidates(I)) .* rangHisto(candidates(I))' ));

for k = 1:length(candidates)
    rangeBin = rang(candidates(I(k)):candidates(I(k))+1);
    idx = find( (rangeBin(1) <= distPV2PC(:)) & (distPV2PC(:) < rangeBin(2)) );
%     for j = 1:sum(~ConnMatPV2PC(idx))
%         if connsNedded(k) > 0
    while (connsNedded(k) > 0) && (sum(~ConnMatPV2PC(idx)))
            tmp = find(~ConnMatPV2PC(idx));
            ConnMatPV2PC(idx(tmp(randperm(length(tmp),1)))) = 1;
            connsNedded(k) = connsNedded(k) - 1 ;
%         end
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

%PV-PV. Connection Probability: O.77 (Gibson, Connors, 1999)
%Eh, o Gibson einai layers 4 and 6, eno o Galarreta, Hestrin, 1999 (same
%issue) einai smoatosensory/visual L5:
% 66% electrical coupling <=80? appart.
% 18% GABAergic connections
% 4.5% Reciprocal GABAergic connections
% ConnMatPV2PV=connectPV2PV(nPV);
[ConnMatPV2PV, gapmat]=connectPV2PV(nPV);
%Label the gap junctions for NEURON (pontiako implementation):
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
AllConnMat = zeros(nAllCells);

% Pyramidals to interneurons:
% PCs to PVs
AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
% PVs connect to all other PVs + autapses 
AllConnMat(pc+1:pv,pc+1:pv) = connectPV2PV(nPV);
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;

% imagesc(AllConnMat)
  

%% Connect pyramidals
% Distance-dependent connection probabilities.
% PC-PC, Non- reciprocal (Nassi)
% BinsPC2PC = 0:30:400;

% Nassis (fer reference)
% total_prob = @(x) 0.22.*exp(-0.006*x);
% recipProbsPC2PC = @(x) 0.15 .* exp(-0.01* x);
% connProbsPC2PC =@(x) ((0.22.*exp(-0.006*x)) - ( 0.12 .* exp(-0.01* x)))*2;
%%
total_prob = @(x) 0.22.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x))); % REMOVE THE MULTIPLICATION FACTOR!!
% SOS remove the multiplication factor because of error in the paper: SOS
upto = 400;
% total_prob = @(x) 0.88.*exp(-0.0052*x);
% recipProbsPC2PC = @(x) 0.48 .* exp(-0.0085* x);
% connProbsPC2PC =@(x) ((0.88.*exp(-0.0052*x)) - (0.48 .* exp(-0.0085* x)))*2;

% overall = imread('overall.jpg');
% reciprocal = imread('reciprocal.jpg');
% nonreciprocal = imread('nonreciprocal.jpg');
fact = 1;%665;
figure;hold on;
% figure;imagesc(flipud(overall));hold on;
% set(gca,'Ydir','normal');
plot(total_prob(0:upto)*fact,'r');

% figure;imagesc(flipud(nonreciprocal));hold on;
% set(gca,'Ydir','normal');
plot(connProbsPC2PC(0:upto)*fact,'k');

% figure;imagesc(flipud(reciprocal));hold on;
% set(gca,'Ydir','normal');
plot(recipProbsPC2PC(0:upto)*fact,'b');
legend({'Overall','Unidirectional','Reciprocal'});

%% Create Networks the Revised way:
myrange = 1:10:max(distPC2PC(:))+10;
tmpDistMat = distPC2PC;
orderedDist = triu(tmpDistMat,1);
orderedDist = orderedDist(logical(triu(ones(nPC,nPC),1)));
distBins = histc( orderedDist, myrange);
% Create Pd:
Pd = zeros(nPC);
for i=1:nPC
    for j=1:nPC
        if (i~=j)
            Pd(i,j)=total_prob(tmpDistMat(i,j));
        else
            Pd(i,j)=0;
        end
    end
end

% Verify complex connectivity:
N = nPC;
E = zeros(N);
start = 2;
for i=1:N
    for j=start:N
        xx=tmpDistMat(i,j);
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
sum(Er(:)) / ((N^2-N)/2)
E = logical(E);

iters = 100;
Es = logical(zeros(N,N,iters));
F_prob = zeros(N,N);
CN_prob = zeros(N,N);
Es(:,:,1) = E;
% connectedPairs_str = [];
% reciprocalPairs_str = [];
connectedPairs_distances_str = {};
reciprocalPairs_distances_str = {};
connectedPairsPercentage_str = [];
reciprocalPairsPercentage_str = [];

p_str_ = zeros(iters,length(myrange));
pr_str_ = zeros(iters,length(myrange));
pnr_str_ = zeros(iters,length(myrange));
GCC = [];
ACC = [];
% fill the first iteration:
t=1;
% connectedPairs_str(t) = (sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1)))) ;
% reciprocalPairs_str(t) = (sum(sum(Es(:,:,t) & Es(:,:,t)'))/2);
connectedPairs_distances_str{t} = tmpDistMat(logical((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1))) ;
reciprocalPairs_distances_str{t} = tmpDistMat(logical((Es(:,:,t) & Es(:,:,t)').*triu(ones(N),1))) ;
connectedPairsPercentage_str(t) = ((sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1)))) / ((N^2 - N)/2));
reciprocalPairsPercentage_str(t) = ((sum(sum(Es(:,:,t) & Es(:,:,t)'))/2) / ((N^2 - N)/2));
orderedConn_str = Es(:,:,t)|Es(:,:,t)';
orderedConn_str = orderedConn_str & logical(triu(ones(size(tmpDistMat)),1));
orderedConn_str = tmpDistMat(logical(orderedConn_str));
p_str = histc( orderedConn_str, myrange)';
p_str_(t,:) = p_str./distBins';
% For reciprocals and non-reciprocals:
pr_str = histc( reciprocalPairs_distances_str{t}, myrange)';
pr_str_(t,:) = pr_str./distBins';
pnr_str = histc( connectedPairs_distances_str{t}, myrange)';
pnr_str_(t,:) = pnr_str./distBins';
% F_prob(:,:,t) = Pd;
% [GCC(t),ACC(t), ~] = clust_coeff(E);

for t=2:iters
    t
    timetaken = tic;
    NCN= m_commonNeighbors(Es(:,:,t-1));
    fprintf('Common neighbors computed in %f seconds.\n',toc(timetaken));

    mean_NCN(t)=mean(NCN(~eye(N)));
    CN_prob=((NCN)./max(NCN(~eye(N))));
%     F_prob(:,:,t)=CN_prob(:,:,t).*Pd;
% To previous update factor ta anevaze sto 8eo. Prepei na min allazoun ta
% possosta twn reciprocals, alla na anadiortanonontai.
    % my damping factor (den doulevei kala..):
%     F_prob(:,:,t) = 0.01 * CN_prob(:,:,t) +  0.99 * Pd;
% Perin's: CN linearly interpolated connection probability values does not
% work either:
% F_prob(:,:,t) = CN_prob(:,:,t) ;
    timetaken = tic;
F_prob = CN_prob .* Pd; %einai simantiko na kanw * Pd !
% F_prob(:,:,t) = F_prob(:,:,t) * (sum(sum(E)) / N^2) * ((N*N)-N)/sum(sum(F_prob(:,:,t)));
F_prob = F_prob * ( (sum(sum(Pd))) / sum(sum(F_prob)) );

% Nassi's solution above is problematic because the individual pair
% probability can be greater than one eliminating completely other pairs p
% in order to keep the sum constant !
    fprintf('F prob updated in %f seconds.\n',toc(timetaken));
    
    timetaken = tic;
    xx=(log(0.22)-log(F_prob))/0.0052;
    xx(xx<0) = 0;
    pr= recipProbsPC2PC(xx);
    pnr=connProbsPC2PC(xx);
    randomNum = rand(nPC);
    tmpRecip = (randomNum<pr) & triu(ones(nPC),1);
    tmpRecip  = tmpRecip | tmpRecip';
    tmpConn = ((randomNum >= pr) & (randomNum < pr+(pnr)));
    tmpConn = tmpConn & triu(ones(nPC),1);
    tmpV = find(tmpConn);
    tmpI = logical(zeros(1,length(tmpV)));
    tmpI(rand(1,length(tmpV))<0.5) = 1;
    tmpConnU = logical(zeros(nPC));
    tmpConnL = logical(zeros(nPC));
    tmpConnU(tmpV(tmpI)) = 1;
    tmpConnL(tmpV(~tmpI)) = 1;
    tmpConnL = tmpConnL';
    Es(:,:,t) = tmpConnU | tmpConnL | tmpRecip;
    fprintf('Es updated in %f seconds.\n',toc(timetaken));

end

for t=2:iters
    t
%     connectedPairs_str(t) = (sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1)))) ;
%     reciprocalPairs_str(t) = (sum(sum(Es(:,:,t) & Es(:,:,t)'))/2);
    connectedPairs_distances_str{t} = tmpDistMat(logical((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1))) ;
    reciprocalPairs_distances_str{t} = tmpDistMat(logical((Es(:,:,t) & Es(:,:,t)').*triu(ones(N),1))) ;
    connectedPairsPercentage_str(t) = ((sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(N),1)))) / ((N^2 - N)/2));
    reciprocalPairsPercentage_str(t) = ((sum(sum(Es(:,:,t) & Es(:,:,t)'))/2) / ((N^2 - N)/2));
    
    
    orderedConn_str = Es(:,:,t)|Es(:,:,t)';
    orderedConn_str = orderedConn_str & logical(triu(ones(size(tmpDistMat)),1));
    orderedConn_str = tmpDistMat(logical(orderedConn_str));
    p_str = histc( orderedConn_str, myrange)';
    p_str_(t,:) = p_str./distBins';
    % For reciprocals and non-reciprocals:
    pr_str = histc( reciprocalPairs_distances_str{t}, myrange)';
    pr_str_(t,:) = pr_str./distBins';
    pnr_str = histc( connectedPairs_distances_str{t}, myrange)';
    pnr_str_(t,:) = pnr_str./distBins';
%     [GCC(t),ACC(t), ~] = clust_coeff(Es(:,:,t));
end

cm=hsv(iters);
cm(1,:) = [0,0,0];
figure;plot(reciprocalPairsPercentage_str ./ connectedPairsPercentage_str);
figure;plot(mean_NCN);
figure;hold on;
for k=iters:-1:1
    plot(myrange,p_str_(k,:),'color',cm(k,:));
end
plot(myrange,total_prob(myrange)*fact,'r','linewidth',3);
figure;hold on;
for k=iters:-1:1
    plot(myrange,pr_str_(k,:),'color',cm(k,:));
end
plot(myrange,recipProbsPC2PC(myrange)*fact,'r','linewidth',3);
figure;hold on;
plot(myrange,p_str_(1,:),'r');
plot(myrange,pr_str_(1,:),'b');
plot(myrange,p_str_(iters,:),'k');
plot(myrange,pr_str_(iters,:),'k');
plot(myrange,recipProbsPC2PC(myrange)*fact,'b','linewidth',3);
plot(myrange,total_prob(myrange)*fact,'r','linewidth',3);
legend({'first','first','last','last','mean','mean'});

prob_conn_rnd_ind = sum(sum(Es(:,:,end))) / N^2;
[Ernd, Ernd_GCC, Ernd_LCC, ~] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind');
sum(sum(Ernd)) / N^2
orderedConn_rnd = Ernd|Ernd';
orderedConn_rnd = orderedConn_rnd & logical(triu(ones(size(tmpDistMat)),1));
sum(orderedConn_rnd(:)) / ((N^2-N)/2)
orderedConn_rnd = tmpDistMat(logical(orderedConn_rnd));
p_rnd = histc( orderedConn_rnd, myrange)';
p_rnd_ = p_rnd./distBins';
% For reciprocals and non-reciprocals:
connectedPairs_distances_rnd = tmpDistMat(logical((Ernd|Ernd').*triu(ones(N),1))) ;
reciprocalPairs_distances_rnd = tmpDistMat(logical((Ernd&Ernd').*triu(ones(N),1))) ;
pr_rnd = histc( reciprocalPairs_distances_rnd, myrange)';
pr_rnd_ = pr_rnd./distBins';
pnr_rnd = histc( connectedPairs_distances_rnd, myrange)';
pnr_rnd_ = pnr_rnd./distBins';
figure;hold on;
plot(myrange,p_rnd_,'b');
plot(myrange,total_prob(myrange),'r','linewidth',3);
figure;hold on;
plot(myrange,pr_rnd_,'color',cm(k,:));
plot(myrange,recipProbsPC2PC(myrange),'r','linewidth',3);
    
sum(sum(Ernd)) / N^2 % random Net
sum(sum(E)) / N^2 % Distance naive Net
sum(sum(Es(:,:,end))) / N^2 % distance Reorganized Net

PC2PC_d = E;
PC2PC_rnd = Ernd;
PC2PC_str = Es(:,:,end);


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
nweights_str = (weights_str/max(weights_str(:)))*2.3;

%% Run Clustering:

mr = -50:1:10; % preference range
tmpidx = cell(1,length(mr));
dpsim = cell(1,length(mr));
[X,Y] = meshgrid(1:N);
similarity_str = zeros((N^2-N)/2,3);
similarity_str(:,1) = X(logical(triu(ones(N),1)));
similarity_str(:,2) = Y(logical(triu(ones(N),1)));
similarity_str(:,3) = CN_str(logical(triu(ones(N),1)));
for k=1:length(mr)
    mr(k)
[tmpidx{k},~,dpsim{k},~]=apclusterSparse(similarity_str,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
[targ_str,~,labels_str] = unique(tmpidx{k}(:,end));

% Random network:
CN_rnd = m_commonNeighbors(PC2PC_rnd);
tmpidx_rnd = cell(1,length(mr));
dpsim_rnd = cell(1,length(mr));
similarity_rnd = zeros((N^2-N)/2,3);
similarity_rnd(:,1) = X(logical(triu(ones(N),1)));
similarity_rnd(:,2) = Y(logical(triu(ones(N),1)));
similarity_rnd(:,3) = CN_rnd(logical(triu(ones(N),1)));
for k=1:length(mr)
    mr(k)
[tmpidx_rnd{k},~,dpsim_rnd{k},~]=apclusterSparse(similarity_rnd,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
[targ_rnd,~,labels_rnd] = unique(tmpidx_rnd{k}(:,end));

figure;plot(cellfun(@(x) x(end), dpsim_rnd))
%% Fill AllConnMat and export it:
% load('C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network\states_07_G.mat')
load('X:\Documents\GitHub\prefrontal-micro\experiment\network\states_07_G.mat')
ID = 12;
SN = 16;
ST = 7;
run = nrun(ID,nPC,SN,ST,1,10000);
pathprefix = 'Z:/data/GliaBackup/BigNetwork_700/';
run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
run.SIMPLIFIED = 1;
% Insert synaptic clustering bias: 0 = clustered, 1 = random (no
% clustering)
run.CLUSTBIAS = 0; % keep synapses in original locations (no rand jittering)
run.FORCECLUSTERING = 0; % Recompute the clusters (overriding the precomputed states)
run.init(states);
run.stimCellsPerCluster = 100; % more cells for better kick and less FF.
run.nPC = nPC;run.nPV=nPV;run.nCB=1;run.nCR=1;run.nAll = run.nPC+run.nPV+run.nCB+run.nCR;
% Create background activity: (keep it simple as in Haider et al., Nature Letter, 2013)
tstop = 20000;
num_exc = 6;
num_inha = 2;
evlen = ((run.nPC*num_exc*3)+(run.nPV*num_inha));
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

% create lognormal weights:
% lognormalW = lognpdf(0:0.01:2.3,0.49,5);
% figure;plot(0:0.01:2.3,lognormalW)
lognormalW = lognrnd(0.49,5,[1,evlen]);

figure;plot(sum(rate))
figure;plot(rate' * lognormalW')

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




stateID = 1;
% run.state_str(:,:) = AllConnMat;
% run.state_str(1:run.nPC,1:run.nPC) = PC2PC_str(:,:,stateID);
% for k = 1: run.nPC
%     run.state_str(k,k,1) = 1;
% end

run.state_str = AllConnMat;
run.state_str(1:nPC, 1:nPC) = PC2PC_str(:,:,stateID);
for k = 1: nPC
    run.state_str(k,k,1) = 1;
end

cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};
    
run.labels_str = labels_str;
run.cellsPerCluster_str = cellsPerCluster_str{stateID};
run.NC_str = length(cellsPerCluster_str{stateID});

% SSS OVERRIDE: stimulate randomly half of the STR network:
run.StimMat_str = [];
[~,iclust] = sort(run.labels_str,'ascend');
for k = 1:7%run.NC_str
    run.StimMat_str(1:100,k) = iclust(100*(k-1)+1:100*k,1);
end

% RND data:
run.state_rnd = AllConnMat;
run.state_rnd(1:nPC,1:nPC) = PC2PC_rnd(:,:,stateID);
% run.state_rnd(1:run.nPC,1:run.nPC) = PC2PC_rnd(:,:,stateID);
for k = 1: nPC
    run.state_rnd(k,k,1) = 1;
end

cellsPerCluster_rnd(stateID) = {histc(labels_rnd,1:size(targ_rnd,1))'};
run.labels_rnd = labels_rnd;
run.cellsPerCluster_rnd = cellsPerCluster_rnd{stateID};
run.NC_rnd = length(cellsPerCluster_rnd{stateID});

% SSS OVERRIDE: stimulate randomly half of the RND network:
run.StimMat_rnd = [];
[~,iclust] = sort(run.labels_rnd,'ascend');
for k = 1:7%run.NC_str
    run.StimMat_str(1:100,k) = iclust(100*(k-1)+1:100*k,1);
end

run.weights_str = (weights_str/max(weights_str(:)))*5.0;
% for random get the mean weight and assign it to every pair:
run.weights_rnd = ones(run.nPC) * mean(run.weights_str(:));

save(sprintf('%s/WorkingLabels_700.mat',pathprefix),'*','-v7.3');
% override: stimulate more cells:
% run.StimMat_str = [];
% for k = 1:run.NC_str
%     tmpstim = 1:75;
%     run.StimMat_str(:,k) = tmpstim(randperm(9));
% end



% export data
% pathprefix = 'C:\Users\steve\Desktop\matlab_files';
run.NC_str=7;run.NC_rnd=7;

run.exportStimulationParameters(pathprefix)
run.exportNetworkParameters(pathprefix)
%initialize some objects first:
run.generateBackgroundActivity();
run.exportBackgroundStimParams(pathprefix)

% export gap connectivity:
% Export STRUCTURED parameter matrices in .hoc file:
fprintf('Exporting importGapParametersSTR.hoc...\n');
fid = fopen([pathprefix,run.SLASH,'importGapParametersSTR.hoc'],'W');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref gapMatrix\n');
fprintf(fid,'gapMatrix = new Matrix(nPVcells, nPVcells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(gapmat)
    for j=1:length(gapmat)
        fprintf(fid,'gapMatrix.x[%d][%d] = %d\n',i-1,j-1, gapmat(i,j));
    end
end
fprintf(fid,'//EOF\n');
fclose(fid);
%%
many_nets=3;

PC2PC_str = zeros(nPC,nPC,many_nets);
PC2PC_d   = zeros(nPC,nPC,many_nets);
PC2PC_rnd = zeros(nPC,nPC,many_nets);

% distPC2PC = ones(nPC,nPC,many_nets);
distPC2PC = repmat(distPC2PC,[1,1,many_nets]);
for k=1:many_nets
    k
    % Change network in space and check what is happening:
%     [PCsomata, distPC2PC(:,:,k)]= CreateRandomNetwork(nPC, cube_dimensions);
    [PC2PC_d(:,:,k), PC2PC_str(:,:,k), coeffs_global(k, 2:3), coeffs_local(k,2:3), ~, pp(k,2:3),Pd]= create_graph_CN(distPC2PC(:,:,k),connProbsPC2PC,recipProbsPC2PC, total_prob);
    prob_conn_rnd_ind =pp(1,3)/2;%0.13; SOS CHANGED THE CONNECTION PROB
    % na to tre3w polles fores k na kratisw afto pou einai pio konta sto
    % input independent probability...
    [PC2PC_rnd(:,:,k), coeffs_global(k, 1), coeffs_local(k,1), pp(k, 1)] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind'); % 1.0=Random graph, 0.0=Watts-Strogatz graph
end

for k=1:many_nets
    k
    connectedPairs_rnd(k) = ((sum(sum((PC2PC_rnd(:,:,k) | PC2PC_rnd(:,:,k)').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2));
    connectedPairs_d(k) = ((sum(sum((PC2PC_d(:,:,k) | PC2PC_d(:,:,k)').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2));
    connectedPairs_str(k) = ((sum(sum((PC2PC_str(:,:,k) | PC2PC_str(:,:,k)').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2));
    reciprocalPairs_rnd(k) = ((sum(sum(PC2PC_rnd(:,:,k) & PC2PC_rnd(:,:,k)'))/2) / ((nPC^2 - nPC)/2));
    reciprocalPairs_d(k) = ((sum(sum(PC2PC_d(:,:,k) & PC2PC_d(:,:,k)'))/2) / ((nPC^2 - nPC)/2));
    reciprocalPairs_str(k) = ((sum(sum(PC2PC_str(:,:,k) & PC2PC_str(:,:,k)'))/2) / ((nPC^2 - nPC)/2));
end



% New MATLAB color values:
% http://www.mathworks.com/help/matlab/graphics_transition/why-are-plot-lines-different-colors.html
y = [mean(connectedPairs_rnd),mean(connectedPairs_d),mean(connectedPairs_str)] * 100;
s = [std([(connectedPairs_rnd);(connectedPairs_d);(connectedPairs_str)]')/sqrt(length(connectedPairs_rnd))]*100;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
% cb = [0    0.4470    0.7410;0.8500    0.3250    0.0980;0.9290    0.6940    0.1250;0 0.6667 0.7373;0.8471 0.4549 0.1765];
cb = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.7773    0.9102    1.0000;
    1.0000    0.8555    0.7969;
    1.0000    0.8867    0.6289];
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', cb(i,:));
    errorbar(i,y(i),s(i),'k');
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'Random', 'Distance', 'Structured'});
title('Connected Pairs (%)')

y = [mean(reciprocalPairs_rnd),mean(reciprocalPairs_d),mean(reciprocalPairs_str)]*100;
s = [std([(reciprocalPairs_rnd);(reciprocalPairs_d);(reciprocalPairs_str)]')/sqrt(length(reciprocalPairs_rnd))]*100;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', cb(i,:));
    errorbar(i,y(i),s(i),'k');
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'Random', 'Distance', 'Structured'});
title('Reciprocal Pairs (%)')
% Find reciprocal percentage from connected pairs.
y = [mean(reciprocalPairs_rnd./connectedPairs_rnd),mean(reciprocalPairs_d./connectedPairs_d),mean(reciprocalPairs_str./connectedPairs_str)] * 100;
s = [std([(reciprocalPairs_rnd./connectedPairs_rnd);(reciprocalPairs_d./connectedPairs_d);(reciprocalPairs_str./connectedPairs_str)]')/sqrt(length(connectedPairs_rnd))]*100;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', cb(i,:));
    errorbar(i,y(i),s(i),'k');
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'Random', 'Distance', 'Structured'});
title('Reciprocal Pairs (%)')



myrange = 1:10:max(distPC2PC(:))+10;
p_bar = [];
p_str_ = [];
p_rnd_ = [];
p_d_ = [];
for k = 1:many_nets
    tmpDistMat = distPC2PC(:,:,k);
    orderedDist = triu(tmpDistMat,1);
    orderedDist = orderedDist(logical(triu(ones(nPC,nPC),1)));
%     orderedDist = orderedDist(orderedDist>0);
    
    orderedConn_str = PC2PC_str(:,:,k)|PC2PC_str(:,:,k)';
    orderedNotConn_str = ~orderedConn_str;
    orderedConn_str = orderedConn_str & logical(triu(ones(size(tmpDistMat)),1));
    orderedNotConn_str = orderedNotConn_str & logical(triu(ones(size(tmpDistMat)),1));
    orderedConn_str = tmpDistMat(logical(orderedConn_str));
    orderedNotConn_str = tmpDistMat(logical(orderedNotConn_str));
    
    orderedConn_rnd = PC2PC_rnd(:,:,k)|PC2PC_rnd(:,:,k)';
    orderedNotConn_rnd = ~orderedConn_rnd;
    orderedConn_rnd = orderedConn_rnd & logical(triu(ones(size(tmpDistMat)),1));
    orderedNotConn_rnd = orderedNotConn_rnd & logical(triu(ones(size(tmpDistMat)),1));
    orderedConn_rnd = tmpDistMat(logical(orderedConn_rnd));
    orderedNotConn_rnd = tmpDistMat(logical(orderedNotConn_rnd));
    
    orderedConn_d = PC2PC_d(:,:,k)|PC2PC_d(:,:,k)';
    orderedNotConn_d = ~orderedConn_d;
    orderedConn_d = orderedConn_d & logical(triu(ones(size(tmpDistMat)),1));
    orderedNotConn_d = orderedNotConn_d & logical(triu(ones(size(tmpDistMat)),1));
    orderedConn_d = tmpDistMat(logical(orderedConn_d));
    orderedNotConn_d = tmpDistMat(logical(orderedNotConn_d));
    
    distBins = histc( orderedDist, myrange);
    % bar(distBins)
    % As a Bernoulli Dist, though p not probability:
    % http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    p_str = histc( orderedConn_str, myrange)';
    p_rnd = histc( orderedConn_rnd, myrange)';
    p_d = histc( orderedConn_d, myrange)';
    q_str = histc( orderedNotConn_str, myrange)';
    q_rnd = histc( orderedNotConn_rnd, myrange)';
    q_d = histc( orderedNotConn_d, myrange)';
    
    % Calculate confidence interval as Agresti-Coull Interval:
    z=1.96;
    n_bar = (numel(tmpDistMat) - size(tmpDistMat,1))/2 + z^2;
    p_str_bar = (p_str+(z^2)/2)/n_bar;
    p_rnd_bar = (p_rnd+(z^2)/2)/n_bar;
    p_d_bar = (p_d+(z^2)/2)/n_bar;
    sem_str = z*sqrt((p_str_bar.*(1-p_str_bar))/n_bar);
    sem_rnd = z*sqrt((p_rnd_bar.*(1-p_rnd_bar))/n_bar);
    sem_d = z*sqrt((p_d_bar.*(1-p_d_bar))/n_bar);
    p_str_(k,:) = p_str./distBins';
    p_rnd_(k,:) = p_rnd./distBins';
    p_d_(k,:) = p_d./distBins';
end

figure;hold on;
for k = 1:many_nets
    plot(1:size(p_str_,2),p_str_(k,:),'color',cb(6,:));
end
for k = 1:many_nets
    plot(1:size(p_rnd_,2),p_rnd_(k,:),'color',cb(4,:));
end
for k = 1:many_nets
    plot(1:size(p_d_,2),p_d_(k,:),'color',cb(5,:));
end
h_str = plot(1:size(p_str_,2),nanmean(p_str_),'Color',cb(3,:),'Linewidth',3);
h_rnd = plot(1:size(p_rnd_,2),nanmean(p_rnd_),'Color',cb(1,:),'Linewidth',3);
h_d = plot(1:size(p_rnd_,2),nanmean(p_d_),'Color',cb(2,:),'Linewidth',3);
h_perin = plot(total_prob(myrange),'k');
set(gca,'XTick',1:5:length(myrange))
set(gca,'XTickLabel',myrange(1:5:end));
legend([h_rnd,h_d,h_str,h_perin],{'Mean Random','Mean Distance', 'Mean Structured','Observed'});
xlabel('Distance (um)');
ylabel('Probability of connection');


% Check weights:
% weights_str = ones(run.nPC);
% weights_rnd = ones(run.nPC);
% for cl = 1:length(cellsPerCluster_str{stateID})
%     idx = find(cl_labels_str(:,stateID) == cl);
%     weights_str(idx,idx) = 1.3; 
% end
% weights_str = weights_str .* PC2PC_str(1:75,1:75,stateID);
% sample random weights from normal distribution:

snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
rang = 0:0.01:8;
m = [0,0,0,0];
alpha = [20,10,5,5];
sigma = [0.45, 0.55, 1.3, 2];
CNsnormal = [];
% Probabilities as in Perin et al.,2011:
% CNsnormal(1,:) = snormal(rang, m(1), alpha(1), sigma(1))*0.012;
% CNsnormal(2,:) = snormal(rang, m(2), alpha(2), sigma(2))*0.014;
% CNsnormal(3,:) = snormal(rang, m(3), alpha(3), sigma(3))*0.02;
% CNsnormal(4,:) = snormal(rang, m(4), alpha(4), sigma(4))*0.02;
% Probabilities changed to replicate mean EPSP amplitude per connections in
% cluster (Perin et al, Fig.6A):
CNsnormal(1,:) = snormal(rang, m(1), alpha(1), sigma(1))*0.012;
CNsnormal(2,:) = snormal(rang, m(2), alpha(2), sigma(2))*0.014;
CNsnormal(3,:) = snormal(rang, m(3), alpha(3), sigma(3))*0.02;
CNsnormal(4,:) = snormal(rang, m(4), alpha(4), sigma(4))*0.02;

% [parmhat,parmci] = lognfit(snormal(rang, m(1), alpha(1), sigma(1))*0.012);
% CNsnormal(1,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(2), alpha(2), sigma(2))*0.014);
% CNsnormal(2,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(3), alpha(3), sigma(3))*0.02);
% CNsnormal(3,:) = lognpdf(rang,parmhat(1),parmhat(2));
% [parmhat,parmci] = lognfit(snormal(rang, m(4), alpha(4), sigma(4))*0.02);
% CNsnormal(4,:) = lognpdf(rang,parmhat(1),parmhat(2));
figure;
plot(CNsnormal(1,:),'b');hold on;
plot(CNsnormal(2,:),'r');hold on;
plot(CNsnormal(3,:),'g');hold on;
plot(CNsnormal(4,:),'k');hold on;
set(gca,'xtick',1:60:length(rang));
set(gca,'xticklabel',rang(1:60:end));

maxCDF = [max(cumsum(CNsnormal(1,:))), max(cumsum(CNsnormal(2,:))), max(cumsum(CNsnormal(3,:))),  max(cumsum(CNsnormal(4,:)))];
weights_str = ones(75,75,many_nets);
for stateID = 1:many_nets
    stateID
    tmp_weights_str = ones(75,75);
    [CN_str] = m_commonNeighbors(PC2PC_str(1:75,1:75,stateID));
    CN_str(CN_str>3) = 3;
    for k=1:4
        sterngthArr = [];
        ctr = 0;
        while ctr < numel(CN_str)%length(find(CN_str==k-1))
            randomSampling = rand(1);
            if randomSampling > maxCDF(k)
                ctr = ctr+1;
                continue;
            end
%             strength = rang(find(histc(rand(1)*max(cumsum(CNsnormal(k,:))),cumsum(CNsnormal(k,:)))));
            strength = rang(find(histc(randomSampling,cumsum(CNsnormal(k,:)))));
            if ~isempty(strength)
                sterngthArr = [sterngthArr, strength(1)];
                ctr = ctr+1;
            end
        end
        tmp_weights_str(CN_str==k-1) = sterngthArr(randperm(length(sterngthArr),length(find(CN_str==k-1))));
    end
    weights_str(:,:,stateID) = tmp_weights_str .* PC2PC_str(1:75,1:75,stateID);
end


averageEPSP = cell(many_nets,30);
for kk = 1:many_nets
    kk
    for k=1:1000
        PCpermut = randperm(run.nPC,6);
        tmpNetId = 1;%ceil(rand(1)*many_nets);
        PCbatch = PC2PC_str(PCpermut,PCpermut,tmpNetId);
        EPSPs = weights_str(PCpermut,PCpermut,tmpNetId);
        EPSPs(EPSPs==0) = [];
        if ~isempty(EPSPs)
            averageEPSP{kk,sum(PCbatch(:))+1} = [averageEPSP{kk,sum(PCbatch(:))+1}, EPSPs];
        end
    end
end
AvrEPSPamplitude = cell(many_nets,1);
for kk = 1:many_nets
    kk
    AvrEPSPamplitude{kk,1} = cell2mat(cellfun(@(x) mean(x),averageEPSP(kk,:),'Uniformoutput',false));
end
EPSPdata = cell2mat(AvrEPSPamplitude);
EPSPstd = nanstd(EPSPdata,1);
EPSPmean = nanmean(EPSPdata,1);
EPSPsem = EPSPstd./sqrt( sum(~isnan(EPSPdata),1) );
figure;errorbar(EPSPmean,EPSPsem);

averageEPSPhist = cellfun(@(x) histc(x,[0:0.1:6]),averageEPSP,'UniformOutput', false);
figure;
for k=1:length(averageEPSPhist)
    subplot(1,length(averageEPSPhist),k);box off;
    barh(averageEPSPhist{k});set(gca,'xticklabel','','yticklabel','');
end
figure;plot(cell2mat(cellfun(@(x) mean(x),averageEPSP,'Uniformoutput',false)))

% Degrees histogam.
degreesHisto_str = [];
degreesHisto_d = [];
degreesHisto_rnd = [];
degreesRange = 1:5:100;
for k =1:many_nets
   degreesHisto_str = [degreesHisto_str ; histcounts( degrees(PC2PC_str(1:nPC,1:nPC,k)),degreesRange)];
   degreesHisto_rnd = [degreesHisto_rnd ; histcounts( degrees(PC2PC_rnd(1:nPC,1:nPC,k)),degreesRange)];
   degreesHisto_d = [degreesHisto_d ; histcounts( degrees(PC2PC_d(1:nPC,1:nPC,k)),degreesRange)];
end
figure;hold on;
plot(degreesHisto_str','Color',cb(6,:));
plot(degreesHisto_rnd','Color',cb(4,:));
plot(degreesHisto_d','Color',cb(5,:));
h_str = plot(mean(degreesHisto_str),'Color',cb(3,:),'LineWidth',3);
h_rnd = plot(mean(degreesHisto_rnd),'Color',cb(1,:),'LineWidth',3);
h_d = plot(mean(degreesHisto_d),'Color',cb(2,:),'LineWidth',3);
% legend([h_str,h_rnd],{'Mean Structured','Mean Random'});
legend([h_rnd,h_d,h_str],{'Mean Random','Mean Distance', 'Mean Structured'});
set(gca,'xtick',1:2:length(degreesRange));
set(gca,'xticklabel',degreesRange(1:2:end));
title('Network Degree');
xlabel('Degree');
ylabel('Frequency');


CN_str = {};
CN_rnd = {};
CN_d = {};
myrange = 1:1:30;
ordNeighHisto_str = [];
ordNeighHisto_rnd = [];
ordNeighHisto_d = [];
for stateID = 1:many_nets
    [CN_str{stateID}] = m_commonNeighbors(PC2PC_str(1:nPC,1:nPC,stateID));
    [CN_rnd{stateID}] = m_commonNeighbors(PC2PC_rnd(1:nPC,1:nPC,stateID));
    [CN_d{stateID}] = m_commonNeighbors(PC2PC_d(1:nPC,1:nPC,stateID));
    orderedNeighbors_str = triu(CN_str{stateID},1);
    orderedNeighbors_rnd = triu(CN_rnd{stateID},1);
    orderedNeighbors_d = triu(CN_d{stateID},1);
    ordNeighHisto_str(stateID,:) = histc( orderedNeighbors_str(logical(triu(ones(size(orderedNeighbors_str)),1))), myrange)';
    ordNeighHisto_rnd(stateID,:) = histc( orderedNeighbors_rnd(logical(triu(ones(size(orderedNeighbors_rnd)),1))), myrange)';
    ordNeighHisto_d(stateID,:) = histc( orderedNeighbors_d(logical(triu(ones(size(orderedNeighbors_d)),1))), myrange)';
end
figure;hold on;
plot(ordNeighHisto_str','Color',cb(6,:));
plot(ordNeighHisto_rnd','Color',cb(4,:));
plot(ordNeighHisto_d','Color',cb(5,:));
h_str = plot(mean(ordNeighHisto_str),'Color',cb(3,:),'LineWidth',3);
h_rnd = plot(mean(ordNeighHisto_rnd),'Color',cb(1,:),'LineWidth',3);
h_d = plot(mean(ordNeighHisto_d),'Color',cb(2,:),'LineWidth',3);
% legend([h_str,h_rnd],{'Mean Structured','Mean Random'});
legend([h_rnd,h_d,h_str],{'Mean Random','Mean Distance', 'Mean Structured'});
set(gca,'xtick',myrange(1:4:end));
set(gca,'xticklabel',myrange(1:4:end));
title('Common Neighbors');
xlabel('Number of Common Neighbors');
ylabel('Frequency');

cellsPerCluster_str = {};
cellsPerCluster_rnd = {};
for stateID = 1:many_nets
    stateID
    % Find nearest neighbors.
    [CN_str] = m_commonNeighbors(PC2PC_str(:,:,stateID));
    [CN_rnd] = m_commonNeighbors(PC2PC_rnd(:,:,stateID));
    
    % keep only the pairs that are connected.
    mergedCN_str = CN_str .* PC2PC_str(:,:,stateID);
    mergedCN_rnd = CN_rnd .* PC2PC_rnd(:,:,stateID);
    
    % Affinity propagation algorithm (Frey, Dueck, 2007):
    % Force different # of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    ClNo = 5;
%     ClNo = -6.3575;
    minClustSize = 9;
    while 1
        labels_str=[];
        % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
        [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
% [idx_str,~,~,~,pref_str]=apcluster(mergedCN_str, ClNo,'dampfact',0.9, ...
%             'convits',200,'maxits',2000,'nonoise');
        [targ_str,~,labels_str] = unique(idx_str);
        if ( min(histc(labels_str,1:ClNo)') >= minClustSize)
            break;
        end
    end
    cl_labels_str(:,stateID)=labels_str;
    [within_str(stateID),between_str(stateID),WB_str(stateID)] = calculateWithinBetween(PC2PC_str(:,:,stateID),labels_str);
%     [within_w_str(stateID),between_w_str(stateID),WB_w_str(stateID)] = calculateWithinBetween(weights_str(:,:,stateID),labels_str);
    cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};
    
    while 1
        labels_rnd=[];
%         NC_rnd=[];
        % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
        [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd, ClNo,0);
% [idx_rnd,~,~,~,pref_rnd]=apcluster(mergedCN_rnd, ClNo,'dampfact',0.9, ...
%             'convits',200,'maxits',2000,'nonoise');
        [targ_rnd,~,labels_rnd] = unique(idx_rnd);
        if (min(histc(labels_rnd,1:ClNo)') > minClustSize)
            break;
        end
    end
    cl_labels_rnd(:,stateID)=labels_rnd;
    [within_rnd(stateID),between_rnd(stateID),WB_rnd(stateID)] = calculateWithinBetween(PC2PC_rnd(:,:,stateID),labels_rnd);
%     [within_w_rnd(stateID),between_w_rnd(stateID),WB_w_rnd(stateID)] = calculateWithinBetween(weights_rnd(:,:,stateID),labels_rnd);
    cellsPerCluster_rnd(stateID)= {histc(labels_rnd,1:size(targ_rnd,1))'};
end

recipPerCluster_str = {};
recipInCluster_str = [];
recipInClusterTotal_str = [];
recipInTotal_str = [];
for k=1:many_nets
    tmpMat = [];
    for l = 1:size(cellsPerCluster_str{k},2)
        lIdx = find(cl_labels_str(:,k)==l);
        tmpMat(l) = sum(sum(PC2PC_str(lIdx,lIdx,k) & PC2PC_str(lIdx,lIdx,k)'))/2;
    end
    recipPerCluster_str{k} = tmpMat;
    recipInTotal_str(k) = sum(sum(PC2PC_str(:,:,k) & PC2PC_str(:,:,k)'))/2;
    recipInCluster_str(k) = sum(recipPerCluster_str{k}) / recipInTotal_str(k);
    recipInClusterTotal_str(k) = sum(recipPerCluster_str{k});
end
recipPerCluster_rnd = {};
recipInCluster_rnd = [];
recipInClusterTotal_rnd = [];
recipInTotal_rnd = [];
for k=1:many_nets
    tmpMat = [];
    for l = 1:size(cellsPerCluster_rnd{k},2)
        lIdx = find(cl_labels_rnd(:,k)==l);
        tmpMat(l) = sum(sum(PC2PC_rnd(lIdx,lIdx,k) & PC2PC_rnd(lIdx,lIdx,k)'))/2;
    end
    recipPerCluster_rnd{k} = tmpMat;
    recipInTotal_rnd(k) = sum(sum(PC2PC_rnd(:,:,k) & PC2PC_rnd(:,:,k)'))/2;
    recipInCluster_rnd(k) = sum(recipPerCluster_rnd{k}) / recipInTotal_rnd(k);
    recipInClusterTotal_rnd(k) = sum(recipPerCluster_rnd{k});
end

y = [mean(recipInCluster_rnd),mean(recipInCluster_str)]*100;
s = [std([(recipInCluster_rnd);(recipInCluster_str)]')/sqrt(length(recipInCluster_str))]*100;
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', cb(i,:));
    errorbar(i,y(i),s(i),'k');
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'Random', 'Structured'});
title('Reciprocals in Clusters (from total)')
ylabel('%');

y = [mean(recipInClusterTotal_rnd),mean(recipInClusterTotal_str)];
s = std([(recipInClusterTotal_rnd);(recipInClusterTotal_str)]')/sqrt(length(recipInClusterTotal_str));
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', cb(i,:));
    errorbar(i,y(i),s(i),'k');
end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'Random', 'Structured'});
title('Total Reciprocals in Clusters (from total)')
ylabel('# of reciprocals');


% xconnections_str = [];
% for stateID = 1:many_nets
%     tmpMat = ones(75);
%     for l = 1:size(cellsPerCluster_str{stateID},2)
%         lIdx = find(cl_labels_str(:,stateID)==l);
%         tmpMat(lIdx,lIdx) = 0;
%     end
%     xconnections_str(stateID) = sum(tmpMat(:));
% end

APrecipInTotal_str = [];
APrecipInCluster_str = {};
APrecipInClusterTotal_str = {};
xconnections_str = [];
labels_str={};
APrecipInTotal_rnd = [];
APrecipInCluster_rnd = {};
APrecipInClusterTotal_rnd = {};
xconnections_rnd = [];
labels_rnd={};
APrecipInTotal_d = [];
APrecipInCluster_d = {};
APrecipInClusterTotal_d = {};
xconnections_d = [];
labels_d={};
aprange = -100:5:30;

for stateID = 1:many_nets
    stateID
    % Find nearest neighbors.
    [CN_str] = m_commonNeighbors(PC2PC_str(:,:,stateID));
    [CN_rnd] = m_commonNeighbors(PC2PC_rnd(:,:,stateID));
    [CN_d] = m_commonNeighbors(PC2PC_d(:,:,stateID));
    
    % keep only the pairs that are connected.
    mergedCN_str = CN_str .* PC2PC_str(:,:,stateID);
    mergedCN_rnd = CN_rnd .* PC2PC_rnd(:,:,stateID);
    mergedCN_d = CN_d .* PC2PC_d(:,:,stateID);
    
    for k = 1:length(aprange)
        [idx_str,~,~,~,pref_str]=apcluster(mergedCN_str, aprange(k),'dampfact',0.9, ...
            'convits',200,'maxits',2000,'nonoise');
        [targ_str,~,labels_str{k,stateID}] = unique(idx_str);
        cellsPerCluster_str = histc(labels_str{k,stateID},1:size(targ_str,1))';
        tmpMat = ones(75);
        tmpMatRecip = [];
        for l = 1:size(cellsPerCluster_str,2)
            lIdx = find(labels_str{k,stateID}==l);
            tmpMat(lIdx,lIdx) = 0;
            tmpMatRecip(l) = sum(sum(PC2PC_str(lIdx,lIdx,stateID) & PC2PC_str(lIdx,lIdx,stateID)'))/2;
        end
        xconnections_str(k,stateID) = sum(tmpMat(:));
        APrecipInTotal_str(k,stateID) = sum(sum(PC2PC_str(:,:,stateID) & PC2PC_str(:,:,stateID)'))/2;
        APrecipInCluster_str{k,stateID} = sum(tmpMatRecip) / APrecipInTotal_str(k);
        APrecipInClusterTotal_str{k,stateID} = sum(tmpMatRecip);
    end
    for k = 1:length(aprange)
        [idx_rnd,~,~,~,pref_rnd]=apcluster(mergedCN_rnd, aprange(k),'dampfact',0.9, ...
            'convits',200,'maxits',2000,'nonoise');
        [targ_rnd,~,labels_rnd{k,stateID}] = unique(idx_rnd);
        cellsPerCluster_rnd = histc(labels_rnd{k,stateID},1:size(targ_rnd,1))';
        tmpMat = ones(75);
        tmpMatRecip = [];
        for l = 1:size(cellsPerCluster_rnd,2)
            lIdx = find(labels_rnd{k,stateID}==l);
            tmpMat(lIdx,lIdx) = 0;
            tmpMatRecip(l) = sum(sum(PC2PC_rnd(lIdx,lIdx,stateID) & PC2PC_rnd(lIdx,lIdx,stateID)'))/2;
        end
        xconnections_rnd(k,stateID) = sum(tmpMat(:));
        APrecipInTotal_rnd(k,stateID) = sum(sum(PC2PC_rnd(:,:,stateID) & PC2PC_rnd(:,:,stateID)'))/2;
        APrecipInCluster_rnd{k,stateID} = sum(tmpMatRecip) / APrecipInTotal_rnd(k);
        APrecipInClusterTotal_rnd{k,stateID} = sum(tmpMatRecip);
    end
    for k = 1:length(aprange)
        [idx_d,~,~,~,pref_d]=apcluster(mergedCN_d, aprange(k),'dampfact',0.9, ...
            'convits',200,'maxits',2000,'nonoise');
        [targ_d,~,labels_d{k,stateID}] = unique(idx_d);
        cellsPerCluster_d = histc(labels_d{k,stateID},1:size(targ_d,1))';
        tmpMat = ones(75);
        tmpMatRecip = [];
        for l = 1:size(cellsPerCluster_d,2)
            lIdx = find(labels_d{k,stateID}==l);
            tmpMat(lIdx,lIdx) = 0;
            tmpMatRecip(l) = sum(sum(PC2PC_d(lIdx,lIdx,stateID) & PC2PC_d(lIdx,lIdx,stateID)'))/2;
        end
        xconnections_d(k,stateID) = sum(tmpMat(:));
        APrecipInTotal_d(k,stateID) = sum(sum(PC2PC_d(:,:,stateID) & PC2PC_d(:,:,stateID)'))/2;
        APrecipInCluster_d{k,stateID} = sum(tmpMatRecip) / APrecipInTotal_d(k);
        APrecipInClusterTotal_d{k,stateID} = sum(tmpMatRecip);
    end
    
end
figure;hold on;
plot(mean(xconnections_str')/5625,'Color',cb(3,:),'Linewidth',2);
plot(mean(cell2mat(APrecipInCluster_str(1:27,1:24))'),'Color',cb(6,:));
plot(mean(xconnections_rnd')/5625,'Color',cb(1,:),'Linewidth',2);
plot(mean(cell2mat(APrecipInCluster_rnd(1:27,1:24))'),'Color',cb(4,:));
plot(mean(xconnections_d')/5625,'Color',cb(2,:),'Linewidth',2);
plot(mean(cell2mat(APrecipInCluster_d(1:27,1:24))'),'Color',cb(5,:));
legend({'Inter-connections Structured','In-cluster Reciprocals Structured','Inter-connections Random','In-cluster Reciprocals Random'});
figure;hold on;
plot(mean(cell2mat(APrecipInClusterTotal_str(1:27,1:24))'),'Color',cb(2,:))
plot(mean(cell2mat(APrecipInClusterTotal_rnd(1:27,1:24))'),'Color',cb(1,:))

figure;hold on;
for stateID = 1:100
    stateID
    for k = 1:length(aprange)
        cellsPerCluster_str = histc(labels_str{k,stateID},1:size(unique(labels_str{k,stateID}),1))';
        scatter(length(cellsPerCluster_str),cell2mat(APrecipInClusterTotal_str(k,stateID)),20,cb(3,:),'o','filled');
        cellsPerCluster_rnd = histc(labels_rnd{k,stateID},1:size(unique(labels_rnd{k,stateID}),1))';
        scatter(length(cellsPerCluster_rnd),cell2mat(APrecipInClusterTotal_rnd(k,stateID)),20,cb(1,:),'o','filled');
        cellsPerCluster_d = histc(labels_d{k,stateID},1:size(unique(labels_d{k,stateID}),1))';
        scatter(length(cellsPerCluster_d),cell2mat(APrecipInClusterTotal_d(k,stateID)),20,cb(2,:),'o','filled');
    end
end

figure;hold on;
for stateID = 1:100
    stateID
    for k = 1:length(aprange)
        cellsPerCluster_str = histc(labels_str{k,stateID},1:size(unique(labels_str{k,stateID}),1))';
        plot(length(cellsPerCluster_str),xconnections_str(k,stateID),'r*');
        cellsPerCluster_rnd = histc(labels_rnd{k,stateID},1:size(unique(labels_rnd{k,stateID}),1))';
        plot(length(cellsPerCluster_rnd),xconnections_rnd(k,stateID),'b*');
    end
end


% Check common neighbour sufficient assumption:
averageConns_str = [];
averageConns_rnd = [];
averageConns_d = [];
for k=1:10000
    PCpermut = randperm(run.nPC,8);
    PCbatch_str = PC2PC_str(PCpermut,PCpermut,ceil(rand(1)*many_nets));
    PCbatch_rnd = PC2PC_rnd(PCpermut,PCpermut,ceil(rand(1)*many_nets));
    PCbatch_d = PC2PC_d(PCpermut,PCpermut,ceil(rand(1)*many_nets));
%     EPSPs = weights_str(PCpermut,PCpermut);
%     EPSPs(EPSPs==0) = [];
%     if ~isempty(PCbatch)
%         averageEPSP{1,sum(PCbatch(:))+1} = [averageEPSP{1,sum(PCbatch(:))+1}, EPSPs];
        averageConns_str = [averageConns_str, sum(PCbatch_str(:))];
        averageConns_rnd = [averageConns_rnd, sum(PCbatch_rnd(:))];
        averageConns_d = [averageConns_d, sum(PCbatch_d(:))];
%     end
end

connsRange = [0:1:64];
figure;hold on;
connsProb_str = histc(averageConns_str,connsRange);
connsProb_str = connsProb_str/k;
plot(connsProb_str,'Color',cb(3,:),'Linewidth',2);
connsProb_rnd = histc(averageConns_rnd,connsRange);
connsProb_rnd = connsProb_rnd/k;
plot(connsProb_rnd,'Color',cb(1,:),'Linewidth',2);
connsProb_d = histc(averageConns_d,connsRange);
connsProb_d = connsProb_d/k;
plot(connsProb_d,'Color',cb(2,:),'Linewidth',2);

% sum(PC2PC_rnd(:))
% sum(PC2PC_d(:))
% sum(PC2PC_str(:))
%No of connected pairs is different because in random case the connections
%get distributed so more pairs are connected.
% fprintf('Random connected pairs are %.3f%%\n',((sum(sum((PC2PC_rnd | PC2PC_rnd').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));
% fprintf('Distance connected pairs are %.3f%%\n',((sum(sum((PC2PC_d | PC2PC_d').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));
% fprintf('Structured connected pairs are %.3f%%\n',((sum(sum((PC2PC_str | PC2PC_str').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2)*100));
% 
% fprintf('Random reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_rnd & PC2PC_rnd'))/2) / ((nPC^2 - nPC)/2)*100));
% fprintf('Distance reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_d & PC2PC_d'))/2) / ((nPC^2 - nPC)/2)*100));
% fprintf('Structured reciprocal pairs are %.3f%%\n',((sum(sum(PC2PC_str & PC2PC_str'))/2) / ((nPC^2 - nPC)/2)*100));

%% Find clusters in networks
WB_str = {};
% WB_str_BC = WB_str;
WB_rnd = {};
% WB_rnd_BC = WB_rnd;
cellsPerCluster_str = {};
for stateID = 1:many_nets
    % Find nearest neighbors.
    [CN_str] = m_commonNeighbors(PC2PC_str(:,:,stateID));
    [ CN_d ] = m_commonNeighbors(PC2PC_d(:,:,stateID));
%     [CN_rnd] = m_commonNeighbors(PC2PC_rnd(:,:,stateID));
    
    % keep only the pairs that are connected.
    mergedCN_str = CN_str .* PC2PC_str(:,:,stateID);
    mergedCN_d   =  CN_d  .* PC2PC_d(:,:,stateID);
%     mergedCN_rnd = CN_rnd .* PC2PC_rnd(:,:,stateID);
    
    % Affinity propagation algorithm (Frey, Dueck, 2007):
    % Force different # of cluster to get as many clusters as you
    % can from the populations with high min cells per cluster:
    ClNo = 5;
    minClustSize = 9;
    Flg = 1;
    while Flg
%         for i=1:5
            labels_str=[];
            NC_str=[];
            % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
            [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
            [targ_str,~,labels_str] = unique(idx_str);
            if ( min(histc(labels_str,1:ClNo)') >= minClustSize)
                Flg = 0;
                ClNo
                break;
            end
%         end
%         ClNo = ClNo - 1;
    end
    cl_labels_str(:,stateID)=labels_str;
    [within_str(stateID),between_str(stateID),WB_str(stateID)] = calculateWithinBetween(PC2PC_str(:,:,stateID),labels_str);
%     [within_w_str(stateID),between_w_str(stateID),WB_w_str(stateID)] = calculateWithinBetween(weights_str(:,:,stateID),labels_str);
    cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};
%     figure('name','Structured')
%     bar(cellsPerCluster_str{:, stateID});
    
%     Flg = 1;
%     while Flg
%         for i=1:20
%             labels_str=[];
%             NC_str=[];
%             % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
%             [idx_d,~,~,~,pref_d]=apclusterK(mergedCN_d, ClNo,0);
%             [targ_d,~,labels_d] = unique(idx_d);
%             if (min(histc(labels_d,1:ClNo)')> minClustSize)
%                 Flg = 0;
%                 ClNo
%                 break;
%             end
%         end
%         ClNo = ClNo - 1
%     end

%     Flg = 1;
%     while Flg
%         for i=1:5
%             labels_rnd=[];
%             NC_rnd=[];
%             % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
%             [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd, ClNo,0);
%             [targ_rnd,~,labels_rnd] = unique(idx_rnd);
%             if (min(histc(labels_rnd,1:ClNo)') > minClustSize)
%                 Flg = 0;
%                 ClNo
%                 break;
%             end
%         end
%         ClNo = ClNo - 1
%     end
% cl_labels_rnd(:,stateID)=labels_rnd;
% [within_rnd(stateID),between_rnd(stateID),WB_rnd(stateID)] = calculateWithinBetween(PC2PC_rnd(:,:,stateID),labels_rnd);
% cellsPerCluster_rnd(stateID)= {histc(labels_rnd,1:size(targ_rnd,1))'};

    
%     cl_labels_d(:,stateID)=labels_d;
%     
%     
%     
%     cellsPerCluster_d(stateID) = {histc(labels_d,1:size(targ_d,1))'};
%     figure('name','Distance')
%     bar(cellsPerCluster_d{:, stateID});
%     
%     cellsPerCluster_rnd(stateID)= {histc(labels_rnd,1:size(targ_rnd,1))'};
%     figure('name','Random')
%     bar(cellsPerCluster_rnd{:, stateID});
    
%     [within_old,between_old,WB_old] = calculateWithinBetween(states.state_str(1:75,1:75,4),states.gparams(4).labels_str)
%     [within_str,between_str,WB_str] = calculateWithinBetween(PC2PC_str,labels_str);
%     [within_d,between_d,WB_d] = calculateWithinBetween(PC2PC_d,labels_d);
%     [within_rnd,between_rnd,WB_rnd] = calculateWithinBetween(PC2PC_rnd,labels_rnd);
%     [within_rnd,within_d,within_str]
%     [between_rnd,between_d,between_str]
%     
%     bar(1,[within_str],'r');hold on;
%     bar(2,[between_str],'k');hold on;
%     legend({'Intra','Inter'})

    
%     temp_str=[];
%     temp_d=[];
%     temp_rnd=[];
%     for i=1:size(targ_str,1)
%         within=PC2PC_str(find(idx_str==targ_str(i)),find(idx_str==targ_str(i)),stateID);
%         between=PC2PC_str(find(idx_str~=targ_str(i)),find(idx_str~=targ_str(i)),stateID);
%         p_str{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_str{stateID,i,2}=sum(sum( between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_str(i)=p_str{stateID,i,1}- p_str{stateID,i,2};
%     end
%     for i=1:size(targ_d,1)
%         within=PC2PC_d(find(idx_d==targ_d(i)),find(idx_d==targ_d(i)),stateID);
%         between=PC2PC_d(find(idx_d~=targ_d(i)),find(idx_d~=targ_d(i)),stateID);
%         p_d{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_d{stateID,i,2}=sum(sum(between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_d(i)=p_d{stateID,i,1}- p_d{stateID,i,2};
%     end
%     for i=1:size(targ_rnd,1)
%         within=PC2PC_rnd(find(idx_rnd==targ_rnd(i)),find(idx_rnd==targ_rnd(i)),stateID);
%         between=PC2PC_rnd(find(idx_rnd~=targ_rnd(i)),find(idx_rnd~=targ_rnd(i)),stateID);
%         p_rnd{stateID,i,1}=sum(sum(within))/((size(within,1)*size(within,1))-size(within,1));
%         p_rnd{stateID,i,2}=sum(sum( between))/((size( between,1)*size( between,1))-size( between,1));
%         temp_rnd(i)=p_rnd{stateID,i,1}- p_rnd{stateID,i,2};
%     end
    
end
% p_dif(stateID,1)= mean(temp_str);
% p_dif(stateID,2)= mean(temp_d);
% p_dif(stateID,3)= mean(temp_rnd);

diags = cellfun(@(x)diag(x),WB_str,'UniformOutput',false);
spread = cellfun(@(x)std(diag(x)),WB_str,'UniformOutput',true);
figure; hold on;
for k = 1:length(diags)
    plot(diags{k});
end

% plot within/between matrices
mat = rand(5);           %# A 5-by-5 matrix of random values from 0 to 1
imagesc(mat);            %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)

textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:5);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

set(gca,'XTick',1:5,...                         %# Change the axes tick marks
        'XTickLabel',{'A','B','C','D','E'},...  %#   and tick labels
        'YTick',1:5,...
        'YTickLabel',{'A','B','C','D','E'},...
        'TickLength',[0 0]);
    

%% Fill AllConnMat and export it:
% load('C:\Users\marouli\Documents\GitHub\prefrontal-micro\experiment\network\states_07_G.mat')
load('X:\Documents\GitHub\prefrontal-micro\experiment\network\states_07_G.mat')
ID = 12;
SN = 16;
ST = 7;
run = nrun(ID,nPC,SN,ST,1,10000);
% pathprefix = 'Z:/data/GliaBackup/FinalClusterRuns/';
% pathprefix = 'Z:/data/GliaBackup/ClusterRunsInhomogenousWeights/';
% pathprefix = 'Z:/data/GliaBackup/';
pathprefix = 'Z:/data/GliaBackup/BigNetwork_700/';
run.pathToHere = 'C:\Users\steve\Documents\GitHub\prefrontal-micro\experiment\network';
run.SIMPLIFIED = 1;
% Insert synaptic clustering bias: 0 = clustered, 1 = random (no
% clustering)
run.CLUSTBIAS = 0; % keep synapses in original locations (no rand jittering)
run.FORCECLUSTERING = 0; % Recompute the clusters (overriding the precomputed states)
run.init(states);
run.stimCellsPerCluster = 100; % more cells for better kick and less FF.

stateID = 1;
% run.state_str(:,:) = AllConnMat;
% run.state_str(1:run.nPC,1:run.nPC) = PC2PC_str(:,:,stateID);
% for k = 1: run.nPC
%     run.state_str(k,k,1) = 1;
% end
run.state_str = AllConnMat;
run.state_str(1:nPC, 1:nPC) = PC2PC_str(:,:,stateID);
for k = 1: nPC
    run.state_str(k,k,1) = 1;
end
% run.state_str = mymat(1:75,1:75);
[CN_str] = m_commonNeighbors(run.state_str(1:nPC,1:nPC));
mergedCN_str = CN_str .* run.state_str(1:nPC,1:nPC);
ClNo = 12;
minClustSize = 9;
Flg = 1;
while Flg
        labels_str=[];
        NC_str=[];
        % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
        [idx_str,~,~,~,pref_str]=apclusterK(mergedCN_str, ClNo,0)
        [targ_str,~,labels_str] = unique(idx_str);
        if ( min(histc(labels_str,1:ClNo)') >= minClustSize)
            Flg = 0;
            ClNo
            break;
        end
end
cl_labels_str(:,stateID)=labels_str;
[within_str(stateID),between_str(stateID),WB_str(stateID)] = calculateWithinBetween(run.state_str(1:nPC,1:nPC),labels_str);
cellsPerCluster_str(stateID) = {histc(labels_str,1:size(targ_str,1))'};
    
run.labels_str = cl_labels_str(:,stateID);
run.cellsPerCluster_str = cellsPerCluster_str{stateID};
run.NC_str = length(cellsPerCluster_str{stateID});
run.StimMat_str = [];
for k = 1:run.NC_str
    tmpstim = find(run.labels_str==k)
    run.StimMat_str(:,k) = tmpstim(randperm(nPC,run.stimCellsPerCluster));
end

% SSS OVERRIDE: stimulate randomly half of the STR network:
run.labels_str = cl_labels_str{25};
% run.cellsPerCluster_str = cellsPerCluster_str{stateID};
run.NC_str = length(unique(cl_labels_str{25}));
run.StimMat_str = [];
[~,iclust] = sort(run.labels_str,'ascend');
for k = 1:7%run.NC_str
    run.StimMat_str(1:100,k) = iclust(100*(k-1)+1:100*k,1);
end
%Dummy RND data:
run.labels_rnd = run.labels_str;
run.NC_rnd = run.NC_str;
run.StimMat_rnd = run.StimMat_str;

run.weights_str = (nweights_str/max(nweights_str(:)))*12;
run.weights_rnd = run.weights_str;
save(sprintf('%s/WorkingLabels_700.mat',pathprefix),'*','-v7.3');
% override: stimulate more cells:
% run.StimMat_str = [];
% for k = 1:run.NC_str
%     tmpstim = 1:75;
%     run.StimMat_str(:,k) = tmpstim(randperm(9));
% end

run.state_rnd = AllConnMat;
run.state_rnd(1:nPC,1:nPC) = PC2PC_rnd(:,:,stateID);
% run.state_rnd(1:run.nPC,1:run.nPC) = PC2PC_rnd(:,:,stateID);
for k = 1: nPC
    run.state_rnd(k,k,1) = 1;
end
[CN_rnd] = m_commonNeighbors(run.state_rnd(1:nPC,1:nPC));
mergedCN_rnd = CN_rnd .* run.state_rnd(1:nPC,1:nPC);

Flg = 1;
while Flg
        labels_rnd=[];
        NC_rnd=[];
        % Nassi, Changed in apclusterK total # of bisections from 20 to 50 (line 75)
        [idx_rnd,~,~,~,pref_rnd]=apclusterK(mergedCN_rnd, ClNo,0)
        [targ_rnd,~,labels_rnd] = unique(idx_rnd);
        if ( min(histc(labels_rnd,1:ClNo)') >= minClustSize)
            Flg = 0;
            ClNo
            break;
        end
end
cl_labels_rnd(:,stateID)=labels_rnd;
[within_rnd(stateID),between_rnd(stateID),WB_rnd(stateID)] = calculateWithinBetween(run.state_rnd(1:nPC,1:nPC),labels_rnd);
cellsPerCluster_rnd(stateID) = {histc(labels_rnd,1:size(targ_rnd,1))'};

run.labels_rnd = cl_labels_rnd(:,stateID);
run.cellsPerCluster_rnd = cellsPerCluster_rnd{stateID};
run.NC_rnd = length(cellsPerCluster_rnd{stateID});
run.StimMat_rnd = [];
for k = 1:run.NC_rnd
%     tmpstim = find(run.labels_rnd==k)
%     run.StimMat_rnd(:,k) = tmpstim(randperm(run.stimCellsPerCluster));
    run.StimMat_rnd(:,k) = randperm(nPC,run.stimCellsPerCluster);
end
% Generate weights matrices:
% Kayaguchi : more potent synapses in reciprocal pairs
% Perin: more potent synapses with more clustering coeff:
% We try to combine them
run.weights_str = ones(nAllCells);
run.weights_rnd = ones(nAllCells);
for cl = 1:run.NC_str
    idx = find(run.labels_str(:,1) == cl);
    run.weights_str(idx,idx) = 1.3; 
end
for cl = 1:run.NC_rnd
    idx = find(run.labels_rnd(:,1) == cl);
    run.weights_rnd(idx,idx) = 1.3;
end

% weights_str = ones(75,75,many_nets);
% stateID = 1
% tmp_weights_str = ones(75,75);
% [CN_str] = m_commonNeighbors(run.state_str(1:75,1:75));
% CN_str(CN_str>maxNeighbors) = maxNeighbors;
% for k=1:maxNeighbors
%     tmp_weights_str(CN_str==k-1) = sterngthArr(k,randperm(maxPrecomps,sum(sum(CN_str==k-1))));
% end
% run.weights_str = tmp_weights_str .* run.state_str(1:75,1:75);
% weightsHisto = run.weights_str;
% weightsHisto(weightsHisto == 0) = [];
% histo = histc(weightsHisto,rang);
% mean(histo)

% export data
pathprefix = 'C:\Users\steve\Desktop\matlab_files';
run.NC_str=7;run.NC_rnd=7;

run.exportStimulationParameters(pathprefix)
run.exportNetworkParameters(pathprefix)

% export gap connectivity:
% Export STRUCTURED parameter matrices in .hoc file:
fprintf('Exporting importGapParametersSTR.hoc...\n');
fid = fopen([pathprefix,run.SLASH,'importGapParametersSTR.hoc'],'W');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref gapMatrix\n');
fprintf(fid,'gapMatrix = new Matrix(nPVcells, nPVcells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(gapmat)
    for j=1:length(gapmat)
        fprintf(fid,'gapMatrix.x[%d][%d] = %d\n',i-1,j-1, gapmat(i,j));
    end
end
fprintf(fid,'//EOF\n');
fclose(fid);

%%
WB={};
NC = [];
CN_str = {};
within = [];
between = [];
searchrange = -100:10:20;
for stateID = 1:many_nets
    [CN_str{stateID}] = m_commonNeighbors(PC2PC_str(1:run.nPC,1:run.nPC,stateID));
        mergedCN_str = CN_str{stateID} .* PC2PC_str(1:run.nPC,1:run.nPC,stateID);
%     for k = 1:length(searchrange)
%         [stateID k]
%         
%         [idx,~,~,~]=apcluster(mergedCN_str,searchrange(k),'dampfact',0.9, ...
%             'convits',200,'maxits',2000,'nonoise');
%         [~,~,labels] = unique(idx);
%         NC(stateID,k) = length(unique(labels));
%         [within(stateID,k),between(stateID,k),WB(stateID,k)] = calculateWithinBetween(PC2PC_str(1:run.nPC,1:run.nPC,stateID),labels);
%     end
end

figure;hold on;
plot(nanmean(within),'c')
plot(nanmean(between),'g')
plot(nanmean(NC)/run.nPC,'color',[0.7,0.7,0.7])
legend({'Within str','Between str','# of Clusters str','Within rnd','Between rnd','# of Clusters rnd'})

myrange = 1:1:30;
p_bar = [];
ordNeighHisto_str = [];
for k = 1:many_nets
    orderedNeighbors = triu(CN_str{k},1);
    ordNeighHisto_str(k,:) = histc( orderedNeighbors(logical(triu(ones(size(orderedNeighbors)),1))), myrange)';
end

%%

WB={};
CN_rnd = {};
NC = [];
within = [];
between = [];
searchrange = -100:10:20;
for stateID = 1:many_nets
    [CN_rnd{stateID}] = m_commonNeighbors(PC2PC_rnd(1:run.nPC,1:run.nPC,stateID));
        mergedCN_rnd = CN_rnd{stateID} .* PC2PC_rnd(1:run.nPC,1:run.nPC,stateID);
%     for k = 1:length(searchrange)
%         [stateID k]
%         
%         [idx,~,~,~]=apcluster(mergedCN_rnd,searchrange(k),'dampfact',0.9, ...
%             'convits',200,'maxits',2000,'nonoise');
%         [~,~,labels] = unique(idx);
%         NC(stateID,k) = length(unique(labels));
%         [within(stateID,k),between(stateID,k),WB(stateID,k)] = calculateWithinBetween(PC2PC_rnd(1:run.nPC,1:run.nPC,stateID),labels);
%     end
end

figure;hold on;
plot(nanmean(within),'r')
plot(nanmean(between),'b')
plot(nanmean(NC)/run.nPC,'k')
legend({'Within','Between','# of Clusters'})

myrange = 1:1:30;
ordNeighHisto_rnd = [];
for k = 1:many_nets
    orderedNeighbors = triu(CN_rnd{k},1);
    ordNeighHisto_rnd(k,:) = histc( orderedNeighbors(logical(triu(ones(size(orderedNeighbors)),1))), myrange)';
end

%%
% Set indices for ease of mind:
pc = size(PCsomata,1);
pv = size(PCsomata,1) + size(PVsomata,1);
   
% Populate the final all-to-all connectivity matrix. In this matrix rows
% are the sources, columns are the targets.
AllConnMat = zeros(nAllCells);
        
% Pyramidals to interneurons:
% PCs to PVs
AllConnMat(1:pc,pc+1:pv) = ConnMatPC2PV;
% PVs connect to all other PVs + autapses 
AllConnMat(pc+1:pv,pc+1:pv) = connectPV2PV(nPV);
% PVs to PCs based on above connectivity (Yuste 2011):
AllConnMat(pc+1:pv,1:pc) = ConnMatPV2PC;
  
% imagesc(AllConnMat)
  
AllWeightMat = ones(size(AllConnMat));
        
% run for the structured network and for the random network as
% well:

RUNS_str = {};
RUNS_rnd = {};

        
        
        %     export network parameters for that experiment:
        cd(mypath);
        system(sprintf('mkdir experiment_%d',aa));
%         exportNetworkPositions(PCsomata,PVsomata,sprintf('%sexperiment_%d/',mypath,aa));
%         % Export pyramidal's clusters ID:
%         exportNetworkCluster([idx_str,idx_rnd],sprintf('%sexperiment_%d/',mypath,aa));
%         
%         save(sprintf('%sexperiment_%d/experiment.mat',mypath,aa),'-v7.3');

     
        Clusters_str(1:nPC,stateID) = cl_labels_str;

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

