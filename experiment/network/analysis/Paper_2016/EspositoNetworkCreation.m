% Built a structured and random network.  
close all; clear all; clc
run.sn = 3;
rng('default') 
rng(run.sn)

% Number of neurons
nPC = 1000;

%  ---- 3D positions ----
inverseDensity =  nthroot(nPC / (75 / 8),3)
cube_dimensions = inverseDensity * 100;
% 3-d points in space. Arguments: #of cells, max seperation distance.
[PCsomata, distPC2PC, distPC2PCwrapped]= CreateRandomNetwork(nPC, cube_dimensions);
distPC2PC = triu(distPC2PC,1) + tril(distPC2PC',-1);
distPC2PCb = distPC2PC;
distPC2PC = triu(distPC2PCwrapped,1) + tril(distPC2PCwrapped',-1);

run.nPC=nPC;

% visualize differences after wrapping (Also do it after clustering?)
tmprn = 0:5:500;
figure;hold on;
title('PC2PC');
plot(tmprn,histc(distPC2PC(logical(triu(ones(nPC),1))),tmprn) / ((nPC^2-nPC)/2),'r');
plot(tmprn,histc(distPC2PCb(logical(triu(ones(nPC),1))),tmprn) / ((nPC^2-nPC)/2),'b');
legend({'No wrapping','Wrapping'});
xlabel('Intersomatic distance (um)');ylabel('Relative frequency');

figure;
scatter3(PCsomata(:,1),PCsomata(:,2),PCsomata(:,3)); 


% Populate the final all-to-all connectivity matrix. In this matrix rows
% are the sources, columns are the targets.
AllConnMat = zeros(run.nPC);

%
% total_prob = @(x) 0.22.*exp(-0.0052*x);
% recipProbsPC2PC = @(x) 0.12 .* exp(-0.0085* x);
% connProbsPC2PC =@(x) ((0.22.*exp(-0.0052*x)) - (0.12 .* exp(-0.0085* x))); % REMOVE THE MULTIPLICATION FACTOR!!
% SOS remove the multiplication factor because of error in the paper: SOS

total_prob = @(x) 0.9.*exp(-0.0052*x);
recipProbsPC2PC = @(x) 0.45 .* exp(-0.0085* x);
connProbsPC2PC =@(x) ((0.9.*exp(-0.0052*x)) - (0.45 .* exp(-0.0085* x))); % REMOVE THE MULTIPLICATION FACTOR!!


upto = cube_dimensions;
fact = 1;%665;
figure;hold on;
plot(total_prob(0:upto)*fact,'r');
plot(connProbsPC2PC(0:upto)*fact,'k');
plot(recipProbsPC2PC(0:upto)*fact,'b');
legend({'Overall','Unidirectional','Reciprocal'});
xlabel('Intersomatic distance (um)');ylabel('Connection probability');

% Create Structured Network:
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


PC2PC_str = Es(:,:,end);


% Average connection probability:
sum(sum( (PC2PC_str|PC2PC_str') & logical(triu(ones(size(PC2PC_str)),1)) )) / ((nPC^2-nPC)/2)

% Get total configurations:
run.configuration_str = AllConnMat;
run.configuration_str(1:run.nPC, 1:run.nPC) = PC2PC_str;
for k = 1: run.nPC
    run.configuration_str(k,k) = 1;
end
% Create weight matrices:

% Extrapolating Perin's findings, more CN, more shift in distribution:
snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
rang = 0:0.01:8;
NCN = m_commonNeighbors(PC2PC_str);
NCN(logical(eye(run.nPC))) = NaN;
maxNeighbors = nanmax(NCN(:));
minNeighbors = nanmin(NCN(:));
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
figure;plot(CNsnormal(minNeighbors:maxNeighbors,:)');hold on;
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

% afto pou 8elw na kratisw einai to lognormal weights distribution
weights_str(:,:) = (tmp_weights_str + rand(size(tmp_weights_str))/100 ).* PC2PC_str;
nweights_str = (weights_str/max(weights_str(:)))*5;
run.weights_str = nweights_str ;

% Umberto methods:

% Calculate local pair-wise relative connection strength:
% can contain nans if weight is zero:
run.Z = abs(run.weights_str - run.weights_str') ./ (run.weights_str + run.weights_str');
% Z upper triangular:
run.Wss = run.Z(find(triu(ones(size(run.Z)),1)));

trianglenum = @(x) (x^2+x)/2
% tmp = find(Wss==0)
% calculate network symmetry:
run.S = 1-(2 * nansum(run.Wss) / (run.nPC*(run.nPC-1)));

myrange =  -0.1:0.01:1.1 ;
figure;bar(myrange(1:end-1), histcounts( run.Wss, myrange ));
ylabel('Count');xlabel('Pairwise relative connection strength (Z)');


%% Save configuration:
% With different serial number
save(fullfile(osDrive(),'Documents','Umberto',sprintf('UmbertoNetworkCreation_SN%d.mat',run.sn)),'run','-v7.3');

connectivityMatrix = PC2PC_str;
weightsMatrix = run.weights_str;
all(find(all([(weightsMatrix(:)==0) , (connectivityMatrix(:)==0)],2)) == find(any([(weightsMatrix(:)==0) , (connectivityMatrix(:)==0)],2)))
