% load% Load data :
Params = importdata('Params.txt');
startAt = 2;
% 
% Params.data(1)=36;
% Params.data(2)=4;
% Params.data(3)=1;
% % check for missing files:
% for c=1:Params.data(1)
%    (2 == exist( ['dendDistances_',sprintf('%d',c),'.txt'], 'file' ) ) ;
%    for pc=1:Params.data(2)
%        (2 == exist( ['soma_',sprintf('%d',pc),'_run_',sprintf('%d',c),'.dat'], 'file' )) ;
%    end
%    for in=1:Params.data(3)
%         (2 == exist( ['insoma_',sprintf('%d',in),'_run_',sprintf('%d',c),'.dat'], 'file' )) ;
%    end
%     
% end

% Load data:
for c=startAt:Params.data(1)
   load( ['dendDistances_',sprintf('%d',c),'.txt'], 'file' )  ;
   for pc=1:Params.data(2)
       load( ['soma_',sprintf('%d',pc),'_run_',sprintf('%d',c),'.txt'], 'file' ) ;
   end
   for in=1:Params.data(3)
        load( ['insoma_',sprintf('%d',in),'_run_',sprintf('%d',c),'.txt'], 'file' ) ;
   end
    
end


MatPersistent=[];
MatNoPersistent=[];
for c=startAt:Params.data(1)
    %     Load distances of current run:
    temp = mean(eval(sprintf('dendDistances_%d',c)),2);
    
    %     Plot every cell's voltage:
    for pc=1:Params.data(2)
        eval(['plot(soma_',sprintf('%d',pc),'_run_',sprintf('%d',c),')']) ;hold on;
        
        allPersistent(pc) = eval(['isPersistent(soma_',sprintf('%d',pc),'_run_',sprintf('%d',c),',2000)']) ;
    end
    
    if all(allPersistent)
        MatPersistent = [MatPersistent,mean(temp)]  ;
        pRuns(c,1) = 1;
        pRuns(c,2) = mean(temp);
    else
        MatNoPersistent = [MatNoPersistent,mean(temp)]  ;
        pRuns(c,1) = 0;
        pRuns(c,2) = mean(temp);
    end
    
    
 pause()
 cla;
end

mean(MatPersistent)
mean(MatNoPersistent)
plot1 = [mean(MatNoPersistent),0;0,mean(MatPersistent)];
bar(plot1);
% [p,h] = signrank(MatPersistent,MatNoPersistent(1:23))
%%
threshold = [40,60];
ranges = 10:20:160;
boxP = zeros(Params.data(2),2);
boxPI = zeros(Params.data(2),length(ranges));
tempP=[];
tempP = find(pRuns(:,1));
tempP = tempP(randperm(length(tempP)));

boxN = zeros(Params.data(2),2);
boxNI = zeros(Params.data(2),length(ranges));
tempN=[];
tempN = find(pRuns(:,1)<1);
tempN = tempN(randperm(length(tempN)));

if(length(tempP)<length(tempN))
    tempN = tempN(1:length(tempP));
else
    tempP = tempP(1:length(tempN));
end

for i=1:length(tempP)
    for pc=1:Params.data(2)
        boxP(pc,1) = boxP(pc,1) + sum(eval( ['(dendDistances_',sprintf('%d',tempP(i)),'(',sprintf('%d',pc),',:)<', sprintf('%d',threshold(1)),')'] )) ;
        boxP(pc,2) = boxP(pc,2) + sum(eval( ['(dendDistances_',sprintf('%d',tempP(i)),'(',sprintf('%d',pc),',:)>', sprintf('%d',threshold(2)),')'] )) ;
        boxPI(pc,:) = boxPI(pc,:) + histc(eval( ['dendDistances_',sprintf('%d',tempP(i)),'(',sprintf('%d',pc),',:)'] ), ranges ) ;
    end
end


for i=1:length(tempN)
    for pc=1:Params.data(2)
        boxN(pc,1) = boxN(pc,1) + sum(eval( ['(dendDistances_',sprintf('%d',tempN(i)),'(',sprintf('%d',pc),',:) <',sprintf('%d',threshold(1)),')'] )) ;
        boxN(pc,2) = boxN(pc,2) + sum(eval( ['(dendDistances_',sprintf('%d',tempN(i)),'(',sprintf('%d',pc),',:) >',sprintf('%d',threshold(2)),')'] )) ;
        boxNI(pc,:) = boxNI(pc,:) + histc(eval( ['dendDistances_',sprintf('%d',tempN(i)),'(',sprintf('%d',pc),',:)'] ), ranges ) ;
    end
end

%plots:
maxAxis = max(max([boxN,boxP]));
figure;bar([boxN(:,1),boxP(:,1)]);axis([-1,Params.data(2)+2,0,maxAxis]);
figure;bar([boxN(:,2),boxP(:,2)]);axis([-1,Params.data(2)+2,0,maxAxis]);

maxAxis = max(max([boxNI,boxPI]));
figure;bar(boxNI');axis([-1,Params.data(2)+2,0,maxAxis]);
figure;bar(boxPI');axis([-1,Params.data(2)+2,0,maxAxis]);
%% Evaluate...
TEMP=zeros(1,8);
VECT=[];
for i=1:Params.data(1)
    for j=1:Params.data(2)
        TEMP = TEMP + histc(eval( ['dendDistances_',sprintf('%d',i),'(',sprintf('%d',j),',:)'] ), ranges ) ;
        VECT = [VECT, eval( ['dendDistances_',sprintf('%d',i),'(',sprintf('%d',j),',:)'] )];
    end
%     bar(TEMP);
%     pause;
%     close all;
%     TEMP=zeros(1,8);
end
bar(TEMP);
%%

cnt=1;
for i=1:100:30000
    temp(cnt) = mean(soma_1_run_1(i:i+99)) ;
    cnt=cnt+1;
end



%%
for i=1:11
    temp = median(eval(sprintf('dendDistances_%d',i)),2); 
    
    name = [' PFC',sprintf('%d',i),];
    figure('NumberTitle','off', 'Name', name);
    expr = ['plot(somaa',sprintf('%d',i),')'];
    subplot(2,2,1); eval(expr);
    title(sprintf('MeanDist=%f',temp(1)));
    
    expr = ['plot(somab',sprintf('%d',i),')'];
    subplot(2,2,2); eval(expr);
    title(sprintf('MeanDist=%f',temp(2)));
    
    expr = ['plot(somac',sprintf('%d',i),')'];
    subplot(2,2,3); eval(expr);
    title(sprintf('MeanDist=%f',temp(3)));
    
    expr = ['plot(somad',sprintf('%d',i),')'];
    subplot(2,2,4); eval(expr);
    title(sprintf('MeanDist=%f',temp(4)));
    
    pause()
    close(gcf);
end
