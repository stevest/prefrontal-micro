close all; clearvars; clc;

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
PC2PC = zeros(nPC,nPC,3);
PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.32); % 1.0=Random graph, 0.0=Watts-Strogatz graph
% PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.2,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
% PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.2,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
PC2PC(:,:,2) = create_graph_DD(distPC2PC,connBinsPC2PC, connProbsPC2PC);
PC2PC(:,:,1) = create_graph_CN(distPC2PC, connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
    recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);

% to random exei parapanw connections!
% calculate Overall Connection Probabilities for each Graph:
OCP=[];
Nr = 10000;
for j=1:3
    cums_N=[];
    cums_U=[];
    cums_R=[];
    cums_O=[];
    ri = ceil(rand(Nr,1)*nPC);
    rj = ceil(rand(Nr,1)*nPC);
    cums_N = find(diag(PC2PC(ri,rj,j)) == 0);
    cums_U = find(xor(diag(PC2PC(ri,rj,j)),diag(PC2PC(rj,ri,j))));
    cums_R = find(diag(PC2PC(ri,rj,j))&diag(PC2PC(rj,ri,j)));
    cums_O = find(diag(PC2PC(ri,rj,j)) == 1);
   
    OCP_N(j) = length(cums_N)/Nr ;
    OCP_U(j) = length(cums_U)/Nr ;
    OCP_R(j) = length(cums_R)/Nr ;
    OCP_O(j) = length(cums_O)/Nr ;
end   

D={};
CC={};
for i=1:3
    D{i} = degrees(PC2PC(:,:,i)) ;
    [~,~,CC{i}] = clust_coeff(PC2PC(:,:,i)) ;
end



% Degrees with increasing prob of connection
figure;hold on;
plot(histc(D{1},1:5:100),'r')  ;
plot(histc(D{2},1:5:100),'g')  ;
plot(histc(D{3},1:5:100),'b')  ;

% Clustering coeff with increasing prob of connection
figure;hold on;
plot(0:0.01:0.1,histc(CC{1},0:0.01:0.1),'r')  ;
plot(0:0.01:0.1,histc(CC{2},0:0.01:0.1),'g')  ;
plot(0:0.01:0.1,histc(CC{3},0:0.01:0.1),'b')  ;
plot([mean(CC{1}),mean(CC{1})],[0,10],'r')  ;
plot([mean(CC{2}),mean(CC{2})],[0,10],'g')  ;
plot([mean(CC{3}),mean(CC{3})],[0,10],'b')  ;

% Number of common neighbors:
[CNI_1 , CNO_1]= commonNeighbors(PC2PC(:,:,1));
[CNI_2 , CNO_2]= commonNeighbors(PC2PC(:,:,2));
[CNI_3 , CNO_3]= commonNeighbors(PC2PC(:,:,3));

figure;hold on;
plot(histc(CNI_1(logical(triu(ones(nPC)))) ,1:20),'r')
plot(histc(CNI_2(logical(triu(ones(nPC)))) ,1:20),'g')
plot(histc(CNI_3(logical(triu(ones(nPC)))) ,1:20),'b')


%%

close all; clearvars; clc;
addpath(genpath('~/Documents/MATLAB/'));
addpath(genpath('~/Documents/GitHub/prefrontal-micro/'));
addpath(genpath('~/Documents/GitHub/neurocommitto/'));
mypath = '/srv/userdata/HOMES/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/';

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
% PC2PC = zeros(nPC,nPC,4);
% PC2PC(:,:,5) = create_graph_WS(distPC2PC,0.2,0.99); % 1.0=Random graph, 0.0=Watts-Strogatz graph
% PC2PC(:,:,4) = create_graph_WS(distPC2PC,0.2,0.5); % 1.0=Random graph, 0.0=Watts-Strogatz graph
% PC2PC(:,:,3) = create_graph_WS(distPC2PC,0.2,0.01); % 1.0=Random graph, 0.0=Watts-Strogatz graph
% PC2PC(:,:,2) = create_graph_DD(distPC2PC,0.2,connBinsPC2PC, connProbsPC2PC);
% PC2PC(:,:,1) = create_graph_CN(distPC2PC,0.2, connBinsPC2PC, connProbsPC2PC, ...%0.07  %0.2 dinei PA
%     recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);

rndGraph = [];
strGraph = [];
ctr = 1;
for i=0.01:0.01:0.9
    i
    strGraph(:,:,ctr) = create_graph_CN(distPC2PC,i, connBinsPC2PC, connProbsPC2PC, ...
        recipBinsPC2PC, recipProbsPC2PC,NconnBins,NincomingProbs,NoutgoingProbs);
    rndGraph(:,:,ctr) = create_graph_WS(distPC2PC,i,0.99);
    ctr = ctr+1;
end

% calculate Overall Connection Probabilities for each Graph:
rndOCP = [];
strOCP = [];
for j=1:90
    cums=[];
    for i=1:1000
        cums(i) = rndGraph(ceil(rand(1)*75),ceil(rand(1)*75),j);
    end
    rndOCP(j) = sum(cums)/length(cums) ;
    
    cums=[];
    for i=1:1000
        cums(i) = strGraph(ceil(rand(1)*75),ceil(rand(1)*75),j);
    end
    strOCP(j) = sum(cums)/length(cums) ;
end
figure;hold on;
plot(rndOCP,'b');
plot(strOCP,'r');
   
   
rndD={};
rndCC={};
for i=1:length(rndGraph)
    i
    rndD{i} = degrees(rndGraph(:,:,i)) ;
    [~,~,rndCC{i}] = clust_coeff(rndGraph(:,:,i)) ;
end

strD={};
strCC={};
for i=1:length(strGraph)
    i
    strD{i} = degrees(strGraph(:,:,i)) ;
    [~,~,strCC{i}] = clust_coeff(strGraph(:,:,i)) ;
end

% D_x_axis = 1:5:100;
% CC_x_axis = 0:0.05:1;
% histD = cellfun(@(x) histc(x,D_x_axis), D, 'uniformoutput',false);
% histCC = cellfun(@(x) histc(x,CC_x_axis), CC, 'uniformoutput',false);
% 
% c_m = cool(length(D));
% 
% figure;hold on;
% for i=1:length(histD)
%     plot(D_x_axis,(histD{i}),'color',c_m(i,:))  ;
% end
% 
% figure;hold on;
% for i=1:length(histCC)
%     plot(CC_x_axis,(histCC{i}),'color',c_m(i,:)) ;
% end


% Degrees with increasing prob of connection
figure;hold on;
for i=1:length(rndD)
    scatter(i,mean(rndD{i}),'b')  ;
    scatter(i,mean(strD{i}),'r')  ;
end

% Clustering coeff with increasing prob of connection
figure;hold on;
for i=1:length(rndCC)
    scatter(i,mean(rndCC{i}),'b') ;
    scatter(i,mean(strCC{i}),'r') ;
end

% Total # of connections with increasing prob of connection
figure;hold on;
for i=1:length(rndCC)
    scatter(i,sum(sum(rndGraph(:,:,i))),'b') ;
    scatter(i,sum(sum(strGraph(:,:,i))),'r') ;
end

% figure;hold on;
% for i=1:length(CC)
%     scatter(i,mean(CC{i})/mean(D{i})) ;
% end
