% Built a structured and random network.  
close all; clear all; clc
aa=1;
rng('default') 
rng(aa)

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
[PVsomata, distPV2PC, distPV2PCwrapped] = CreateCubeNetworkPV(cube_dimensions, nPV, PCsomata);
distPV2PCb = distPV2PC;
distPV2PC = distPV2PCwrapped;

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

%% Create Networks the Revised way:
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
iters = 1000;
Es = false(nPC,nPC,iters);
F_prob = zeros(nPC,nPC);
CN_prob = zeros(nPC,nPC);
Es(:,:,1) = E;


sumPd = sum(sum(Pd));
for t=2:iters
    t
    NCN= m_commonNeighbors(Es(:,:,t-1));
    
%     mean_NCN(t)=mean(NCN(~eye(nPC)));
    CN_prob=((NCN)./max(NCN(~eye(nPC))));
    
    F_prob = CN_prob .* Pd; %einai simantiko na kanw * Pd !
    F_prob = F_prob * ( sumPd / sum(F_prob(:)) );
    
    timetaken = tic;
    xx=(log(0.22)-log(F_prob))/0.0052;
    xx(xx<0) = 0;
    pr= recipProbsPC2PC(xx);
    pnr=connProbsPC2PC(xx);
    randomNum = rand(nPC);
    % na dw min exw kanei kamia malakia...
    % reciprocals einai swsta:
    tmpRecip = (randomNum<pr) & triu(ones(nPC),1);
    tmpRecip  = tmpRecip | tmpRecip';
    % unidir einai swsta:
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
    % ola kala ws edw:
    Es(:,:,t) = tmpConnU | tmpConnL | tmpRecip;
    
    fprintf('Es updated in %f seconds.\n',toc(timetaken));
    
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
%%
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
% connectedPairs_str(t) = (sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1)))) ;
% reciprocalPairs_str(t) = (sum(sum(Es(:,:,t) & Es(:,:,t)'))/2);
connectedPairs_distances_str{t} = distPC2PC(logical((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1))) ;
reciprocalPairs_distances_str{t} = distPC2PC(logical((Es(:,:,t) & Es(:,:,t)').*triu(ones(nPC),1))) ;
connectedPairsPercentage_str(t) = ((sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2));
reciprocalPairsPercentage_str(t) = ((sum(sum(Es(:,:,t) & Es(:,:,t)'))/2) / ((nPC^2 - nPC)/2));
orderedConn_str = Es(:,:,t)|Es(:,:,t)';
orderedConn_str = orderedConn_str & logical(triu(ones(size(distPC2PC)),1));
orderedConn_str = distPC2PC(logical(orderedConn_str));
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
%     connectedPairs_str(t) = (sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1)))) ;
%     reciprocalPairs_str(t) = (sum(sum(Es(:,:,t) & Es(:,:,t)'))/2);
    connectedPairs_distances_str{t} = distPC2PC(logical((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1))) ;
    reciprocalPairs_distances_str{t} = distPC2PC(logical((Es(:,:,t) & Es(:,:,t)').*triu(ones(nPC),1))) ;
    connectedPairsPercentage_str(t) = ((sum(sum((Es(:,:,t) | Es(:,:,t)').*triu(ones(nPC),1)))) / ((nPC^2 - nPC)/2));
    reciprocalPairsPercentage_str(t) = ((sum(sum(Es(:,:,t) & Es(:,:,t)'))/2) / ((nPC^2 - nPC)/2));
    
    
    orderedConn_str = Es(:,:,t)|Es(:,:,t)';
    orderedConn_str = orderedConn_str & logical(triu(ones(size(distPC2PC)),1));
    orderedConn_str = distPC2PC(logical(orderedConn_str));
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

prob_conn_rnd_ind = sum(sum(Es(:,:,end))) / nPC^2;
[Ernd, Ernd_GCC, Ernd_LCC, ~] = create_graph_WS(nPC,prob_conn_rnd_ind,'ind');
sum(sum(Ernd)) / nPC^2
orderedConn_rnd = Ernd|Ernd';
orderedConn_rnd = orderedConn_rnd & logical(triu(ones(size(distPC2PC)),1));
sum(orderedConn_rnd(:)) / ((nPC^2-nPC)/2)
orderedConn_rnd = distPC2PC(logical(orderedConn_rnd));
p_rnd = histc( orderedConn_rnd, myrange)';
p_rnd_ = p_rnd./distBins';
% For reciprocals and non-reciprocals:
connectedPairs_distances_rnd = distPC2PC(logical((Ernd|Ernd').*triu(ones(nPC),1))) ;
reciprocalPairs_distances_rnd = distPC2PC(logical((Ernd&Ernd').*triu(ones(nPC),1))) ;
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
    
sum(sum(Ernd)) / nPC^2 % random Net
sum(sum(E)) / nPC^2 % Distance naive Net
sum(sum(Es(:,:,end))) / nPC^2 % distance Reorganized Net

PC2PC_d = E;
PC2PC_rnd = Ernd;
PC2PC_str = Es(:,:,end);
