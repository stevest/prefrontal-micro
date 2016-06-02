%% Compute relative frequency of network states:
close all;clear all;clc;
% load(fullfile(osDrive(),'Documents','Glia','NetworkCreation_SN4.mat'));
load(fullfile(osDrive(),'Documents','Glia','NetworkCreation_SN2.mat'));
%Which cluster is stimulated in each configuration:
stc_rnd = 2;
stc_str = 4;
VARPID = 75;
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('SAVENMDAupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_PID%d_spikes.mat',stc_rnd-1, run.sn,VARPID)));
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('SAVENMDAupdatedStimGABAb01NEWBGST_Ss10c%d_SN%d_PID%d_spikes.mat',stc_str-1,run.sn,VARPID)));
run.tstop = 10000;
run.nruns = 100;

% List of stimulated/non-stimulated cells in each configuration:
sc_rnd = run.stimulatedCells_rnd{stc_rnd};
nsc_rnd = find(~ismember(1:700,run.stimulatedCells_rnd{stc_rnd}));
sc_str = run.stimulatedCells_str{stc_str};
nsc_str = find(~ismember(1:700,run.stimulatedCells_str{stc_str}));


%% Array of windows Q to apply:
Qseq = [100,50,40,30,20,10,8,6,4];

% Na e3etazw ono to network, mono to stimulated k mono to recruited 'H ola
% ta ypolloipa kyttara (ola ane3artita apo to an einai recruited).

% Only stimulated cluster:

% Qseq = [100,50,40,30,20,10,8,6,4];
configuration = 'str';
eval( sprintf('st = batch_%s_spikes(sc_%s,:);',configuration,configuration) );
% eval( sprintf('st = batch_%s_spikes;',configuration) );
eval( sprintf('stc = stc_%s;', configuration) );
[~, ~, ~] = createVoteState(run, Qseq, st, stc-1, configuration, 'save');
% Qseq = [100,50,40,30,20,10,8,6,4];
configuration = 'rnd';
eval( sprintf('st = batch_%s_spikes(sc_%s,:);',configuration,configuration) );
eval( sprintf('stc = stc_%s;', configuration) );
[~, ~, ~] = createVoteState(run, Qseq, st, stc-1, configuration, 'save');

%% Get most active cells in each configuration.
FFr_rnd = zeros(run.nPC,run.nruns);
FFr_str = zeros(run.nPC,run.nruns);

for ru=1:size(batch_rnd_spikes,2)
    for c = 1:run.nPC
        st_rnd = batch_rnd_spikes{c,ru};
        st_rnd = st_rnd(st_rnd>1500);
        FFr_rnd(c,ru) = length(st_rnd)/((run.tstop-1500)/1000);
    end
end
for ru=1:size(batch_str_spikes,2)
	for c=1:run.nPC
        st_str = batch_str_spikes{c,ru};
        st_str = st_str(st_str>1500);
        FFr_str(c,ru) = length(st_str)/((run.tstop-1500)/1000);
    end
end

% Sort by mean activity (per cell):
[MEAN_rnd,IDX_rnd] = sort(mean(FFr_rnd,2),'descend');
[MEAN_str,IDX_str] = sort(mean(FFr_str,2),'descend');

cmmax = max(max([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,:)]));
cmmin = min(min([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,:)]));

figure;imagesc(FFr_rnd(IDX_rnd,:));
caxis manual
caxis([cmmin cmmax]);
cm = hot(1000);
cm(1,:) = [0,0,0];
colormap(cm);title('Random');
figure;imagesc(FFr_str(IDX_str,:));
caxis manual
caxis([cmmin cmmax]);
colormap(cm);title('Structured');



%% Compare states across configurations:
VARPID = 75;
close all;
nProminentStatesCheck = 10;
cm = lines(nProminentStatesCheck);
rlowessWidth = 40;
uptocoactive = 10; % plot up to 10 coactive cells.
figure(5);hold on;
for Q = 50
    configuration = 'rnd';
    eval( sprintf('stc = stc_%s', configuration) );
    load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('PID%d_iNMDA_Qanalysis_stimulatedClusterOnly_GABAb01_SN2',VARPID),...
        sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Q)));
    delayRange = ceil(1500/Q):run.tstop/Q ;
    prominentStates = mean(voteState(:,delayRange),2);
    [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
    nActivePC_rnd = cellfun(@(x) length(regexp(x,'1','match')), U );

    frequencyPerActive = {}; gPC={};
    for k=1:uptocoactive 
        tmp = find(nActivePC_rnd == k);
        gPC{k,1} = ones(length(tmp),1)*k;
        frequencyPerActive{k,1} = prominentStates(tmp);
    end

    % get a better visualization for these points:
    figure(2);scatter(cell2mat(gPC),cell2mat(frequencyPerActive),'.');hold on;
    % Show points only for state traces:
    figure(2);scatter(nActivePC_rnd(maxfreqidx(2:nProminentStatesCheck)),prominentStates(maxfreqidx(2:nProminentStatesCheck)),'o');
    
    ylim2_rnd = get(gca,'YLim');
    title(sprintf('%s Q=%d (n(U)=%d)',configuration,Q,size(voteState,1)));
    ylabel('Relative Frequency (%)');xlabel('Coactive neurons');
    figure(5);plot(maxfreqstates(1:nProminentStatesCheck));
    
    title(sprintf('Q=%d',Q));
    ylabel('Relative Frequency (%)');xlabel('States ID (sorted)');
    figure(1);hold on;
    for k=2:nProminentStatesCheck
        plot(smooth(voteState(maxfreqidx(k),delayRange)',rlowessWidth,'rlowess'),'linewidth',2);
    end
    ylim_rnd = get(gca,'YLim');
    title(sprintf('%s Q=%d (n(U)=%d)',configuration,Q,size(voteState,1)));
    ylabel('Relative Frequency (%)');xlabel('Time (in Q windows)');
    
    % plot state intra weight:
    stateMeanW_rnd = [];
    eval( sprintf('stimulatedCells = sc_%s;',configuration) );
    for k=2:nProminentStatesCheck
        [~, S] = regexp(U{maxfreqidx(k)},'1','match');
        stateCells = stimulatedCells(S)';
        if length(stateCells) > 1
            eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
            stateMeanW_rnd(k) = mean(tmpW(~eye(length(stateCells))));
        end
    end
    AllStateMeanW_rnd = {};
    for k=2:uptocoactive
        tmp = find(nActivePC_rnd == k);
        for kk = 1:length(tmp)
            [~, S] = regexp(U{tmp(kk)},'1','match');
            stateCells = stimulatedCells(S)';
            if length(stateCells) > 1
                eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
                AllStateMeanW_rnd{k,kk} = mean(tmpW(~eye(length(stateCells))));
            end
        end
    end
    for k=2:uptocoactive
        tmp2 = cell2mat(AllStateMeanW_rnd(k,:));
        figure(6);scatter(ones(1,length(tmp2))*k,tmp2,'.');hold on;
        figure(6);scatter(k,mean(tmp2),'r+');
    end
    title(sprintf('%s Q=%d',configuration,Q));
    ylabel('Synaptic Weight');xlabel('Coactive neurons');
    
%     % plot state intra syn location:
%     PSLoc_rnd=[];
%     medianLocation_rnd = cellfun(@median,loccummat_rnd);
%     AllStateMeanLoc_rnd = cell(1,uptocoactive);
%     for k=2:uptocoactive
%         tmp = find(nActivePC_rnd == k);
%         for kk = 1:length(tmp)
%             [~, S] = regexp(U{tmp(kk)},'1','match');
%             PSLocTMP = find(ismember(maxfreqidx(2:nProminentStatesCheck),tmp(kk)))+1; % we exclude 1 coactive 
%             if PSLocTMP
%                 PSLoc_rnd = [PSLoc_rnd ; [k,kk,PSLocTMP]];
%             end
%             stateCells = stimulatedCells(S)';
%             if length(stateCells) > 1
%                 tmpL = medianLocation_rnd(stateCells,stateCells);
%                 AllStateMeanLoc_rnd{k} = [AllStateMeanLoc_rnd{k} nanmedian(tmpL(~eye(length(stateCells))))];
%             end
%         end
%     end
%     for k=2:uptocoactive
%         tmp2 = AllStateMeanLoc_rnd{k};
%         figure(8);scatter(ones(1,length(tmp2))*k,tmp2,'k.');hold on;
%         figure(8);scatter(k,nanmedian(tmp2),'r+', 'linewidth',2);
%     end
%     for k=1:size(PSLoc_rnd,1)
%         figure(8);scatter(PSLoc_rnd(k,1),AllStateMeanLoc_rnd{PSLoc_rnd(k,1)}(PSLoc_rnd(k,2)),'ro', 'linewidth',2);
%     end
%     title(sprintf('%s Q=%d',configuration,Q));
%     ylabel('Incomming synapse location (normalized)');xlabel('Coactive neurons');
    
    
    configuration = 'str';
    eval( sprintf('stc = stc_%s', configuration) );
    load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('PID%d_iNMDA_Qanalysis_stimulatedClusterOnly_GABAb01_SN2',VARPID),...
        sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Q)));
    delayRange = ceil(1500/Q):run.tstop/Q ;
    prominentStates = mean(voteState(:,delayRange),2);
    [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
    nActivePC_str = cellfun(@(x) length(regexp(x,'1','match')), U );
    
    frequencyPerActive = {}; gPC={};
    for k=1:uptocoactive
        tmp = find(nActivePC_str == k);
        gPC{k,1} = ones(length(tmp),1)*k;
        frequencyPerActive{k,1} = prominentStates(tmp);
    end

    % get a better visualization for these points:
    figure(4);scatter(cell2mat(gPC),cell2mat(frequencyPerActive),'.');hold on;
    % Show points only for state traces:
    figure(4);scatter(nActivePC_str(maxfreqidx(2:nProminentStatesCheck)),prominentStates(maxfreqidx(2:nProminentStatesCheck)),'o');
    
    ylim2_str = get(gca,'YLim');
    title(sprintf('%s Q=%d (n(U)=%d)',configuration,Q,size(voteState,1)));
    ylabel('Relative Frequency (%)');xlabel('Coactive neurons');
    figure(5);plot(maxfreqstates(1:nProminentStatesCheck));
    
    title(sprintf('Q=%d',Q));
    ylabel('Relative Frequency (%)');xlabel('States ID (sorted)');
    figure(3);hold on;
    for k=2:nProminentStatesCheck
        plot(smooth(voteState(maxfreqidx(k),delayRange)',rlowessWidth,'rlowess'),'linewidth',2);
    end
    ylim_str = get(gca,'YLim');
    title(sprintf('%s Q=%d (n(U)=%d)',configuration,Q,size(voteState,1)));
    ylabel('Relative Frequency (%)');xlabel('Time (in Q windows)');
    
    % plot state intra weight:
    stateMeanW_str = [];
    eval( sprintf('stimulatedCells = sc_%s;',configuration) );
    for k=2:nProminentStatesCheck
        [~, S] = regexp(U{maxfreqidx(k)},'1','match');
        stateCells = stimulatedCells(S)';
        if length(stateCells) > 1
            eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
            stateMeanW_str(k) = mean(tmpW(~eye(length(stateCells))));
        end
    end
    AllStateMeanW_str = {};
    for k=2:uptocoactive
        tmp = find(nActivePC_str == k);
        for kk = 1:length(tmp)
            [~, S] = regexp(U{tmp(kk)},'1','match');
            stateCells = stimulatedCells(S)';
            if length(stateCells) > 1
                eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
                AllStateMeanW_str{k,kk} = mean(tmpW(~eye(length(stateCells))));
            end
        end
    end
    for k=2:uptocoactive
        tmp2 = cell2mat(AllStateMeanW_str(k,:));
        figure(7);scatter(ones(1,length(tmp2))*k,tmp2,'.');hold on;
        figure(7);scatter(k,mean(tmp2),'r+');
    end
    title(sprintf('%s Q=%d',configuration,Q));
    ylabel('Synaptic Weight');xlabel('Coactive neurons');
    
%     % plot state intra syn location:
%     PSLoc_str=[];
%     medianLocation_str = cellfun(@median,loccummat_str);
%     AllStateMeanLoc_str = cell(1,uptocoactive);
%     for k=2:uptocoactive
%         tmp = find(nActivePC_str == k);
%         for kk = 1:length(tmp)
%             [~, S] = regexp(U{tmp(kk)},'1','match');
%             PSLocTMP = find(ismember(maxfreqidx(2:nProminentStatesCheck),tmp(kk)))+1; % we exclude 1 coactive 
%             if PSLocTMP
%                 PSLoc_str = [PSLoc_str ; [k,kk,PSLocTMP]];
% %                 for kkk=1:length(PSLocTMP)
% %                     PSLoc_str{k} = {PSLoc_str(k) {k kk}};
% %                 end
%             end
%             stateCells = stimulatedCells(S)';
%             if length(stateCells) > 1
%                 tmpL = medianLocation_str(stateCells,stateCells);
% %                 AllStateMeanLoc_str{k,kk} = nanmean(tmpL(~eye(length(stateCells))));
%                 AllStateMeanLoc_str{k} = [AllStateMeanLoc_str{k} nanmedian(tmpL(~eye(length(stateCells))))];
%             end
%         end
%     end
%     for k=2:uptocoactive
%         tmp2 = AllStateMeanLoc_str{k};
%         figure(9);scatter(ones(1,length(tmp2))*k,tmp2,'k.');hold on;
%         figure(9);scatter(k,nanmedian(tmp2),'r+', 'linewidth',2);
%     end
%     for k=1:size(PSLoc_str,1)
%         figure(9);scatter(PSLoc_str(k,1),AllStateMeanLoc_str{PSLoc_str(k,1)}(PSLoc_str(k,2)),'ro', 'linewidth',2);
%     end
%     title(sprintf('%s Q=%d',configuration,Q));
%     ylabel('Incomming synapse location (normalized)');xlabel('Coactive neurons');

    
    figure(5);legend({'Rnd','Str'});
    ylim_max = max([ylim_str(2),ylim_rnd(2)]);
    figure(1);set(gca,'YLim',[0, ylim_max]);
    figure(3);set(gca,'YLim',[0, ylim_max]);
    ylim2_max = max([ylim2_str(2),ylim2_rnd(2)]);
    figure(2);set(gca,'YLim',[0, ylim2_max]);
    figure(4);set(gca,'YLim',[0, ylim2_max]);
    
%     % Plot active cluster's activity in each configuration:
%     batch_rnd_spikes
%     m = size(st,2) ; %No of chains (m)
%     n = run.tstop ; %No of itterations (n)
%     N = size(st,1);
%     Qr = floor(n / Q) ; % length of reshaped spiketrain array
%     wst = reshape(st,1,[])';
%     winst = cell(m*N,1);
%     parfor c = 1:m*N
%         if ~isempty(wst{c})
%             tmpst = zeros(1,n);
%             tmpwinst = zeros(1,Qr);
%             tmpst(round(wst{c})) = 1;
%             for k=1:Qr
%                 tmpwinst(k) = any( tmpst( ((((k)-1)*Q)+1):((k) * Q) ) ,2) ;
%             end
%            winst{c} =  tmpwinst;
%         else
%             winst{c} = zeros(1,Qr);
%         end
%     end
%     wspktrain =  cell2mat( cellfun(@(x) reshape(x,1,1,[]),reshape(winst,N,[]),'uniformoutput',false) );
%     clear wst winst;
% 
%     figure(6);hold on;
%     for kk=1:length(sc_rnd)
%         spikes = batch_rnd_spikes{sc_rnd(kk),20};
%         scatter(spikes,ones(length(spikes),1)*kk);
%     end
%     figure(7);hold on;
%     for kk=1:length(sc_str)
%         spikes = batch_str_spikes{sc_str(kk),20};
%         scatter(spikes,ones(length(spikes),1)*kk);
%     end
end

%% Check location of synapses across runs:
loccummat_rnd = cell(700,700);
for k=1:run.nruns
    k
    clear synapticDelays synapticLocations 
    load(sprintf('synapticLocDel_Rs10c%d_SN%d_r%d.mat',stc_rnd-1,run.sn,k-1));
    for ii=1:700
        for jj=1:700
            loccummat_rnd{ii,jj} = [loccummat_rnd{ii,jj} synapticLocations{ii,jj}];
        end
    end
%     cumMat(:,:,k) = cellfun(@mean,synapticLocations);
end
figure;imagesc(cellfun(@median,loccummat_rnd))

loccummat_str = cell(700,700);
for k=1:run.nruns
    k
    clear synapticDelays synapticLocations 
    load(sprintf('synapticLocDel_Ss10c%d_SN%d_r%d.mat',stc_str-1,run.sn,k-1));
    for ii=1:700
        for jj=1:700
            loccummat_str{ii,jj} = [loccummat_str{ii,jj} synapticLocations{ii,jj}];
        end
    end
end
figure;imagesc(cellfun(@median,loccummat_str))
    
%%
% VARPID = 25;
% stc_rnd = 2;
% inmdaSum = zeros(700,700,50);
% for ru=1:50
%     load(fullfile(osDrive(),'Documents','Glia',sprintf('SAVENMDAupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_inmda_r%d_PID%d.mat',stc_rnd-1, run.sn,ru-1,VARPID)));
% 
%     blah = cellfun(@length, inmdaBatch);
%     inmdaIDX = find(blah);
% 
%     for k = inmdaIDX'
%         [y,x] = ind2sub(size(inmdaBatch),k);
%         inmda = inmdaBatch{y,x}(1500:end);
%         inmda(inmda>0)=0;
%         inmdaSum(y,x,ru)=-sum(inmda);
%     end
% end

% figure;imagesc(std(inmdaSum(sc_rnd,sc_rnd,:),0,3))
% figure;imagesc(sum(inmdaSum(sc_rnd,sc_rnd,:),3))

inmdaSTD_R_PID25 = std(inmdaSum_R_PID25(sc_rnd,sc_rnd,:),0,3);
inmdaSUM_R_PID25 = sum(inmdaSum_R_PID25(sc_rnd,sc_rnd,:),3);

inmdaSTD_S_PID25 = std(inmdaSum_S_PID25(sc_str,sc_str,:),0,3);
inmdaSUM_S_PID25 = sum(inmdaSum_S_PID25(sc_str,sc_str,:),3);

figure;hold on;
scatter(inmdaSTD_R_PID25(:),inmdaSUM_R_PID25(:),'.k');
title('RND');
ylabel('SUM');xlabel('STD');
    
    
figure;hold on;
scatter(inmdaSTD_S_PID25(:),inmdaSUM_S_PID25(:),'.k');
title('STR');
ylabel('SUM');xlabel('STD');
    


figure;
for ru=1:49
    imagesc(inmdaSum(sc_rnd,sc_rnd,ru));
    pause();
end

%% Parse above data:
% Na apofasisw ti na kanw me afti tin analysi: Ti apo ola afta pou exw
% kanei toso kairo a3izei na graftei k na to exw, pou exei ginei poutana o
% kwdikas.
configuration = 'str';
nProminentStatesCheck = 10;
for Q = Qseq
    eval( sprintf('stc = stc_%s', configuration) );
    load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab','Qanalysis_stimulatedClusterOnly_GABAb01_SN3',sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Q)));
    %Get most frequent states by smooth data:
    delayRange = ceil(1500/Q):run.tstop/Q ;
    prominentStates = zeros(size(voteState,1),1);
    for kk=1:size(voteState,1)
        prominentStates(kk) = mean(voteState(kk,delayRange)); % SSS changed!
    end
    [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
%     largeFreq = maxfreqstates(maxfreqstates~=0);
    %plot 
    figure;plot(maxfreqstates(1:nProminentStatesCheck));
    title(sprintf('Q=%d',Q));
    ylabel('Relative Frequency (%)');xlabel('States ID (sorted)');
    figure;hold on;
    for k=2:nProminentStatesCheck
        plot(smoothed_states(maxfreqidx(k),delayRange))
    end
    title(sprintf('Q=%d',Q));
    ylabel('Relative Frequency (%)');xlabel('Time (in Q windows)');
end



    
%%

% n=run.tstop;
% Qr = floor(n / Q)
% maxstates = max(smoothed_states,[],2);
% [maxfreqstates,maxfreqidx] = sort(maxstates,'descend') ;
% largeFreq = maxfreqstates(maxfreqstates~=0);
% idx = maxfreqidx(1:length(largeFreq));
% % No of subplot columns:
% ncols=7;
% for s=1:1%ceil(length(largeFreq)/(nstates))
%     % Den to pistevw oti grafw ton parakatw kwdika:
%     if (nstates + nstates*(s-1)) <=length(largeFreq)
%         ftmpind = (1:nstates) + nstates*(s-1);
%     else
%         ftmpind =  floor(length(largeFreq)/nstates)*nstates+1:length(largeFreq) ;
%     end
%     if ncols>length(ftmpind)
%         ncols = length(ftmpind);
%     end
%     % Get Absolute index to states array U:
%     mystates = idx(ftmpind);
%     % leave outside the zero state when calculating maxy:
%     maxy = ones(length(mystates),1)*max(maxfreqstates(~(maxfreqidx==1)));
%     maxy(maxfreqidx==1) = maxfreqstates(maxfreqidx==1);
%     
%     nrows = ceil(length(mystates)/ncols);
%     figure('Name',sprintf('%s, Cluster: %d, Q: %d, From: %d TH:%f',upper(configuration),stc,Q,nstates*(s-1),0),'Position', [ 63           1        1858        1003]);
%     for k=1:length(mystates)
%         
%         subplot(nrows,ncols,k);hold on;
%         plot(voteState(mystates(k),:),'x');
%         plot(smoothed_states(mystates(k),:),'r', 'linewidth',2);% SSS changed to moving for speed!
%         % Set axis limits:
%         xlim([0,Qr]);
%         ylim([0,maxy(k)]);
%         title(sprintf('Active State Neurons: %d',length(regexp(U{mystates(k)},'1','match')) ));
%         if ~(mod(k-1,ncols)==0)
%             set(gca, 'YTick', []);
%         end
%         if ~( (((nrows-1)*ncols)+1)<=k )
%             set(gca, 'XTick', []);
%         end
%         
%     end
% end



% for ii=1:100
% % Plot instantaneous firing frequency:
% IFRArray = zeros(nPC,run.tstop);
% for k=1:nPC
%     spikes = st{k,ii};
%     if ~isempty(spikes)
%     IFRArray(k,:) = instantfr(spikes,linspace(0,run.tstop,run.tstop));
%     end
% end
% 
% cm = jet(round(max(IFRArray(:)))*100);
% cm(1,:)=[0,0,0];
% figure;imagesc(IFRArray(run.StimMat_str(:),:));
% colormap(cm);
% pause();
% close;
% end
%% Yiota: distribution of states for STR vs RND (across different Qs):
tstop = 20000;
% run.nruns = 100;
m = run.nruns ; %No of chains (m)
n = tstop; %No of itterations (n)

Qseq = 100;%[2:1:50,55:5:100];

U={};
% voteState = cell(length(Qseq),7);
for Qs = 1%:length(Qseq)
    Q = Qseq(Qs) ; % simple window (ms)
    Qr = floor(n / Q) ; % length of wspiketrain
    n_over_b = 20 ; % No of batches
    b = floor(Qr/n_over_b); % Batches' length (b)
    if b < 2
        error('Batches'' length must be greater than one sample!');
    end
    
    s = (((1:Qr)-1) * Q) + 1 ;
    e = ((1:Qr) * Q) ;
    
    for stc=1%:NC
        % Calculate spike trains ONLY for the above selected cells:
        tic;
        wspktrain = zeros(run.nPC,m,Qr);
%         for ru = 1:m
%             spktrain = zeros(run.nPC,n);
%             for c=1:run.nPC % must be row vector!!
%                 if ~isempty(RUNS{1,stc}{c,ru})
%                     spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
%                 end
%             end
%             for k=1:Qr
%                 wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
%             end
%         end
        for ru = 1:m
            spktrain = zeros(run.nPC,n);
            for c=1:run.nPC % must be row vector!!
                if ~isempty(st_str{c,ru})
                    spktrain(c,round(st_str{c,ru}')) = 1;
                end
            end
            for k=1:Qr
                wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
            end
        end
        
        fprintf('Generating spiketrain took: %fs\n',toc);
        %Work with string representations:
        wstates = reshape(wspktrain,700,[])';
        wstrstates = cell(length(wstates),1);
        for kk=1:length(wstates)
            wstrstates{kk}=sprintf('%i',wstates(kk,:));
        end
        U{stc,Qs} = unique(wstrstates);
        voteState = zeros(size(U{stc,Qs},1),Qr);
        for qr=1:Qr
%             B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
            B = cell(100,1);
            tmp = squeeze(wspktrain(:,:,qr))';
            for kk=1:100
                B{kk}=sprintf('%i',tmp(kk,:));
            end
            [~,locStates] = ismember(B,U{stc,Qs});
            if any(locStates==0)
                error('This should never happen..');
            end
            voteState(:,qr) = accumarray([1:size(U{stc,Qs},1), locStates']',[voteState(:,qr)', ones(1,length(locStates))]) ;
        end
        
        
        %             Work with binary representations:
        tic;F = sort(bi2de(reshape(wspktrain,700,[])'));
        fprintf('Generating sorted binary representations took: %fs\n',toc);
        
        tic;U{stc,Qs} = unique(F); % briskei ligotera values apo tin allh me8odo WTF
        % alla me to parakatw testing ola fainetai na leitourgoun kanonika:
        % test binary transformations: for k=1:700;dearray(k) = bi2de(de2bi(k));end;all(diff(dearray)==1)
        % epishs to unique fainetai na bgainei to idio, opws kai na to kanw
%         TMP = {};
%         UTMP = {};
%         for qr=1:Qr
%             TMP{qr} = bi2de(squeeze(wspktrain(:,:,qr))');
%             UTMP{qr} = unique(TMP{qr});
%         end
%         UUTMP = unique(cell2mat(UTMP'));
%         all( sort(UUTMP) == sort(G{stc}) )
        % Epomenws to most prominent state gia ena cluster mporei na min
        % anhkei entos tou cluster???

        fprintf('Generating unique states took: %fs\n',toc);
        tic;
%         voteState{stc,Qs} = zeros(size(U{stc,Qs},1),Qr);
        voteState = zeros(size(U{stc,Qs},1),Qr);
        for qr=1:Qr
            B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
            [~,locStates] = ismember(B,U{stc,Qs});
            voteState(:,qr) = accumarray([1:size(U{stc,Qs},1), locStates']',[voteState(:,qr)', ones(1,length(locStates))]) ;
        end
        %             voteState{stc,Qs} = voteState{stc,Qs}./m;
        fprintf('Votting table took: %fs\n',toc);
        
%         % Na sortarw ta uniqe states analoga me to Hamming dist tous me to
%         % stim state.kai na dw posa kyttara exei to ka8ena.
%         tmpArray = zeros(length(U{stc,Qs}),4);
%         tmpstimvec = zeros(nPC,1);
%         tmpstimvec(run.StimMat_str(:,stc)) = 1;
%         for k=1:length(U{stc,Qs})
%             tmpArray(k,1) = pdist2(padarray(de2bi(U{stc,Qs}(k)),[0,700-length(de2bi(U{stc,Qs}(k)))],'post'),tmpstimvec','hamming');
%             tmpArray(k,2) = sum(de2bi(U{stc,Qs}(k)));
%             tmpArray(k,3) = max(voteState{stc,Qs}(k,:));
%             tmpArray(k,4) = mean(voteState{stc,Qs}(k,:));
%         end
%         only=1:1500;
%         % Sort by Hamming dist:
%         [~,idx1]=sort(tmpArray(:,1),'ascend');
%         %Sort by # of active cells:
%         [~,idx2]=sort(tmpArray(:,2),'descend');
%         figure;plot(tmpArray(idx1(only),1));
%         figure;plot(tmpArray(idx1(only),2));
%         figure;plot(tmpArray(idx1(only),3));
%         figure;plot(tmpArray(idx1(only),4));
%         %Parapanw fainetai na exoume polla states me >5% freq kai >4 cells
        
        % Ta 3anakanw edw, mipws ta exei ftysei h matlab:
        %TO GET DETAILED STATE:
        stm = ceil(1500/Q);
        [whymax,idx] = sort(max(voteState(:,stm:end),[],2),'descend');
%         tmps = size(voteState{stc,Qs},2)-stm;
%         [whymax,idx] = sort(mean(voteState{stc,Qs}(:,floor(tmps /10)*9:end),2),'descend');
        data = voteState(idx,stm:end)./m;
        tic;
        maxstates = zeros(size(voteState,1),1);
        fprintf('Generating max of smoothed states...\n');
        for kk=1:size(voteState,1)
            tmp = smooth(voteState(kk,:)',100,'rlowess');
            maxstates(kk) = max(tmp(floor(Qr/2):Qr));
        end
        fprintf('Generating max of smoothed states took: %fs\n',toc);
        [maxfreqstates,maxfreqidx] = sort(maxstates,'descend') ;
        largeFreq = maxfreqstates(maxfreqstates~=0);
        idx = maxfreqidx(1:length(largeFreq));
%         %Test proof data:
%         cm = jet(max(max(data))*100);
%         cm(1,:) = [0,0,0];
%         figure;imagesc(sort(sort(data,1,'descend'),2,'descend'));
%         colormap(cm);
        
        maxdata = max(data,[],2);
%         maxdata = mean(data(:,floor(tmps /10)*9:end),2);

        es = 29;%length(find(maxdata>=0.1))+1;
        ss = 2;
        if ss:es
            mystates = ss:es;
        else 
            mystates = ss:ss+1;
        end

%         figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Q),'Position', [100, 100, 640, 480]);hold on;plot(maxdata)
%         title(sprintf('Q: %d',Q));
%         if ss < es
%             axis([ss es+1 0 maxdata(ss)])
%         else
%             axis([ss ss+1 0 maxdata(ss)])
%         end

%         if length(mystates)>3
%         figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Q),'Position', [740, 100, 1049, 895]);
%         else
%         figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Q),'Position', [740, 100, 1049, 357]);
%         end
%         for k=1:length(mystates)
%             if length(mystates)<3
%                 tria = length(mystates) ;
%             else
%                 tria = 3;
%             end
%             subplot(ceil(length(mystates)/tria),tria,k);hold on;
%             plot(data(mystates(k),:));
% %             gaussFilter = gausswin(Qr/10);
% %             gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
% 
% %             tmp = conv(data(mystates(k),:), gaussFilter,'same');
%             tmp = smooth(data(mystates(k),:)',Qr/10,'rlowess');
% %             tmp = tmp(1:end-length(gaussFilter));
%             plot(tmp,'r', 'linewidth',2);
%             title(sprintf('Active State Neurons: %d',sum(de2bi(U{stc,Qs}(idx(k)))) ));
%         end
        
        figure('Name',sprintf('RND, Cluster: %d, Q: %d',stc,Q));
        for k=1:length(mystates)
            tria=7;
            subplot(ceil(length(mystates)/tria),tria,k);hold on;
            plot(data(mystates(k),:),'x');
            tmp = smooth(data(mystates(k),:)',100,'rlowess');
            plot(tmp,'r', 'linewidth',2);
            title(sprintf('Active State Neurons: %d',sum(de2bi(U{stc,Qs}(idx(k)))) ));
        end

        
%         pause(0.5);
%         imageData = screencapture(0,[0,0,1920,1080]);
%         imwrite(imageData,sprintf('RND_Cluster%d_Q%d.jpg',stc,Q),'jpg','Quality',100);
%         close all;
%         
    end
end
clear spktrain wspktrain;
% save('voteState_str.mat','voteState','-v7.3');
%% plot interpolated data from above:
% sort with sum of all frequencies (for all time points):
[tmp,idx] = sort(sum(voteState{1,1}/m,2),'descend');
plot(tmp(2:end));

cm = hsv(max(voteState{1,1}(:)));
cm(1,:) = [0,0,0];
imagesc(voteState{1,1}(idx,:));
colormap(cm)

% generate interpolated table to surf later:
interpCell = [];
for stc = 1:7
    for Qseq = 1:length(2:1:50)
        interpCell(1:500,Qseq,stc) = interp1(linspace(0,1,length(U{stc,Qseq})),U{stc,Qseq},linspace(0,1,500))';
    end
end

surf( interpCell(1:end-1,:,1)/max(max(interpCell(:,:,1))),'linestyle','none');

figure;plot(U{1,1}(1:end-1));
figure;plot(U{1,49}(1:end-1),'r');
%%
m = run.nruns ; %No of chains (m)
n = tstop; %No of itterations (n)

% Qseq = [2,5,10,15,20,30,40,80,150,200,500];
% Qseq = [2,3,4,5,6,7];
Qseq = 2:2:50;
Nseq = 1;

% voteState = cell(length(Qseq),7);
for Qs = 1:length(Qseq)
    for Ns = 1:length(Nseq)
        Q = Qseq(Qs) ; % simple window (ms)
        Qr = floor(n / Q) ; % length of wspiketrain
        n_over_b = 20 ; % No of batches
        b = floor(Qr/n_over_b); % Batches' length (b)
        if b < 2
            error('Batches'' length must be greater than one sample!');
        end
        
        s = (((1:Qr)-1) * Q) + 1 ;
        e = ((1:Qr) * Q) ;
        
%         clusterStatesAll = cell(1,7);
        for stc=1:NC
            % Calculate spike trains ONLY for the above selected cells:
            wspktrain = zeros(run.nPC,m,Qr);
            for ru = 1:m
                spktrain = zeros(run.nPC,n);
                for c=1:run.nPC % must be row vector!!
                    if ~isempty(RUNS_str{1,stc}{c,ru})
                        spktrain(c,round(RUNS_str{1,stc}{c,ru}.spikes')) = 1;
                    end
                end
                for k=1:Qr
                    wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
                end
            end
            
            % Panos: check for frequency of states across runs:
            %             Work with binary representations:
            tic;F = sort(bi2de(reshape(wspktrain,700,[])'));toc;
            tic;G = unique(F);toc; % briskei ligotera values apo tin allh me8odo WTF
            
            voteState{Qs,stc} = zeros(size(G,1),Qr);
            for qr=1:Qr
                qr
                B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
                [~,locStates] = ismember(B,G);
                % Gamhse...
                voteState{Qs,stc}(:,qr) = accumarray([1:size(G,1), locStates']',[voteState{Qs,stc}(:,qr)', ones(1,length(locStates))]) ;
            end
            voteState{Qs,stc} = voteState{Qs,stc}./m;
            
            % %             % Break to batches to avoid out of memory:
            % %             %chunks = 1:700*run.nruns:numel(wspktrain);
            % %             chunks = 1:size(wspktrain,3);
            % %             clusterStates = {};
            % %             for chunk = 1:length(chunks)
            % %                 [clusterStates{chunk}, ~, ~] = unique(squeeze(wspktrain(:,:,chunk))', 'rows');
            % %             end
            % %             clusterStatesUnion = [];
            % %             clusterStatesUnion = cell2mat(cellfun(@(x)x',clusterStates,'uniformoutput',false));
            % %             [clusterStatesAll{stc}, ~, ~] = unique(clusterStatesUnion', 'rows');
            % %             clearvars clusterStatesUnion clusterStates
            % %             % Panos: check for frequency of states across runs:
            % %             voteState{Qs,stc} = zeros(size(clusterStatesAll{stc},1),Qr);
            % %             for qr=1:Qr
            % %                 qr
            % %                 % the clusterStatesAll must be wrong!
            % %                     %[~,locStates] = ismember(squeeze(wspktrain(:,qr,:))', clusterStatesAll{stc}, 'rows');
            % %                     [~,locStates] = ismember(squeeze(wspktrain(:,:,qr))', clusterStatesAll{stc}, 'rows');
            % %                     % Gamhse...
            % %                     voteState{Qs,stc}(:,qr) = accumarray([1:size(clusterStatesAll{stc},1), locStates']',[voteState{Qs,stc}(:,qr)', ones(1,length(locStates))]) ;
            % %             end
            % %             voteState{Qs,stc} = voteState{Qs,stc}./m;
        end
        
        
    end
end
% save('X:\Documents\Glia\str_states_distributions_7clusters.mat','voteState_str','-v7.3')
% save('X:\Documents\Glia\rnd_states_distributions_7clusters.mat','voteState_rnd','-v7.3')
% Construct blurring window.
windowWidth = int16(50);
halfWidth = windowWidth / 2;
gaussFilter = gausswin(200);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

%% Plot distribution of states across time:
[~,sidx] = sort(sum(voteState_str{20,1},2),'descend');
surf(voteState_str{20,1}(sidx(2:100),:),'linestyle','none')
%%
% %             posa diaforetika states gia ka8e qr? DONE
% %             a3onas apo 0-1 DONE
% %             meta to stimulation DONE
% %             windows apo 2-10
% %             smothing DONE
% %             Panos: posa states emfanizontai (max) megalyterh syxnothta across runs
% Ena tropo na kanw visualize apo 2-50ms Q posa einai ta prominent states
% (se frequency prominent) kai gi aftes to smooth subplot. posa einai sta
% states. h sygklish tou ka8enos. diaferei rnd VS str. Spesificity?
% oti ta most prominent a8roizoun sto 100 konta
% na balw to 00 state
% na kanw visualize thn diafora twn states (specificity hamming distance)

% load detailed states builded above:
% voteState_str = voteState_str';
% voteState_str = cellfun(@(x) x./m,voteState_str,'uniformoutput',false);

Qseq = 2:1:50;
for Qs = 1:49
    
    Q = Qseq(Qs) ; % simple window (ms)
    Qr = floor(n / Q) ; % length of wspiketrain
    n_over_b = 20 ; % No of batches
    b = floor(Qr/n_over_b); % Batches' length (b)
    if b < 2
        error('Batches'' length must be greater than one sample!');
    end
    s = (((1:Qr)-1) * Q) + 1 ;
    e = ((1:Qr) * Q) ;
    
    for stc = 1%:NC
        wspktrain = zeros(run.nPC,m,Qr);
        for ru = 1:m
            spktrain = zeros(run.nPC,n);
            for c=1:run.nPC % must be row vector!!
                if ~isempty(RUNS{1,stc}{c,ru})
                    spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
                end
            end
            for k=1:Qr
                wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
            end
        end
        clear spktrain;
        tic;F = sort(bi2de(reshape(wspktrain,700,[])'));toc; % na bebaiw8w oti einai swstos o tropos aftos (EINAI)
        tic;G{stc} = unique(F);toc; % briskei ligotera values apo tin allh me8odo WTF (SWSTO)
        tmpvoteState = zeros(size(G{stc},1),Qr);
        for qr=1:Qr
            B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
            [~,locStates] = ismember(B,G{stc});
            tmpvoteState(:,qr) = accumarray([1:size(G{stc},1), locStates']',[tmpvoteState(:,qr)', ones(1,length(locStates))]) ;
        end
   
        %How many neuros are active per state:
        statesActiveNeurons = zeros(length(G{stc}),1);
        for k=1:length(G{stc})
            statesActiveNeurons(k) = sum(de2bi(G{stc}(k)));
        end
        
        %TO GET DETAILED STATE:
        stm = ceil(1500/Qseq(Qs));
        [whymax,idx] = sort(max(tmpvoteState(:,stm:end),[],2),'descend');

%         % How is the no of neurons correlated with frequency of occurance:
%         figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Qseq(Qs)),'Position', [287   287   560   420]);
%         scatter(statesActiveNeurons(2:end)',max(tmpvoteState(2:end,stm:end)/m,[],2));title('Max Freq');ylabel('Relative Frequency');xlabel('# of active neurons per state');
%         figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Qseq(Qs)),'Position', [849   287   560   420]);
%         scatter(statesActiveNeurons(2:end)',mean(tmpvoteState(2:end,stm:end)/m,2));title('Mean Freq');ylabel('Relative Frequency');xlabel('# of active neurons per state');

        
% %         cm = hsv(max(max(tmpvoteState(:,stm:end))));
%         cm = hsv(2);
%         cm(1,:) = [0,0,0];
%         figure;imagesc(wspktrain(:,:,1000))
% %         figure;imagesc(tmpvoteState(:,stm:end));
%     %     figure;imagesc(data);
%         colormap(cm)


    %ORIGINAL:            
    %     % calculate stimulation overhead:
    %     stm = ceil(1500/Qseq(Qs));
    %     
    %     [whymax,idx] = sort(max(voteState_str{Qs,stc}(:,stm:end),[],2),'descend');
    %     data = voteState_str{Qs,stc}(idx,stm:end)./m;
        data = tmpvoteState(idx,stm:end)./m;

        % figure;plot(sum(logical(data)))
        maxdata = max(data,[],2);
        es = length(find(maxdata>=0.1))+1;
        ss = 2;
        if ss:es
            mystates = ss:es;
        else 
            mystates = ss:ss+1;
        end

        figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Qseq(Qs)),'Position', [100, 100, 640, 480]);hold on;plot(maxdata)
        title(sprintf('Q: %d',Qseq(Qs)));
        if ss < es
            axis([ss es+1 0 maxdata(ss)])
        else
            axis([ss ss+1 0 maxdata(ss)])
        end

        if length(mystates)>3
        figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Qseq(Qs)),'Position', [740, 100, 1049, 895]);
        else
        figure('Name',sprintf('STR, Cluster: %d, Q: %d',stc,Qseq(Qs)),'Position', [740, 100, 1049, 357]);
        end
        for k=1:length(mystates)
            if length(mystates)<3
                tria = length(mystates) ;
            else
                tria = 3;
            end
            subplot(ceil(length(mystates)/tria),tria,k);hold on;
            plot(data(mystates(k),:));
            gaussFilter = gausswin(Qr/10);
            gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

            tmp = conv(data(mystates(k),:), gaussFilter,'same');
            tmp = tmp(1:end-length(gaussFilter));
            plot(tmp,'r', 'linewidth',2);
            title(sprintf('Active State Neurons: %d',sum(de2bi(G{stc}(idx(k)))) ));

            %plot(conv(data(mystates(k),:), gaussFilter,'same'),'r', 'linewidth',2);
    %         plot(smoothts(data(mystates(k),:),'g',500,5000))
        end
        
        
        pause(0.5);
        imageData = screencapture(0,[0,0,1920,1080]);
%         imwrite(imageData,sprintf('STR_Cluster%d_Q%d.jpg',stc,Qseq(Qs)),'jpg','Quality',100);
        imwrite(imageData,sprintf('STR_FreqVSNeurons_Cluster%d_Q%d.jpg',stc,Qseq(Qs)),'jpg','Quality',100);
        close all;
        
    end % for stc
end

%% Calculate Hamming distance between clusters:

m = 97%run.nruns ; %No of chains (m)
n = tstop; %No of itterations (n)

Qseq = 5;
Nseq = 1;

for Qs = 1:length(Qseq)
%     for Ns = 1:length(Nseq)
        Q = Qseq(Qs) ; % simple window (ms)
        Qr = floor(n / Q) ; % length of wspiketrain
        n_over_b = 20 ; % No of batches
        b = floor(Qr/n_over_b); % Batches' length (b)
        if b < 2
            error('Batches'' length must be greater than one sample!');
        end
        s = (((1:Qr)-1) * Q) + 1 ;
        e = ((1:Qr) * Q) ;

        BI = {};
        DE = {};
        G = {};
        mystates = {};
        %Panos: mono ena cluster gia arxh:
        for stc = 1:NC
            wspktrain = zeros(run.nPC,m,Qr);
            for ru = 1:m
                spktrain = zeros(run.nPC,n);
                for c=1:run.nPC % must be row vector!!
                    if ~isempty(RUNS_rnd{1,stc}{c,ru})
                        spktrain(c,round(RUNS_rnd{1,stc}{c,ru}.spikes')) = 1;
                    end
                end
                for k=1:Qr
                    wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
                end
            end
            clear spktrain;
            tic;F = sort(bi2de(reshape(wspktrain,700,[])'));toc; % na bebaiw8w oti einai swstos o tropos aftos
            tic;G{stc} = unique(F);toc; % briskei ligotera values apo tin allh me8odo WTF
            
%             tic;U{stc,Qs} = unique(F); % briskei ligotera values apo tin allh me8odo WTF
            tmpvoteState = zeros(size(G{stc},1),Qr);
            for qr=1:Qr
                B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
                [~,locStates] = ismember(B,G{stc});
                tmpvoteState(:,qr) = accumarray([1:size(G{stc},1), locStates']',[tmpvoteState(:,qr)', ones(1,length(locStates))]) ;
            end
        
        
            stm = ceil(1500/Qseq(Qs));
%             [whymax,idx] = sort(max(voteState_str{Qs,stc}(:,stm:end),[],2),'descend');
            [whymax,idx] = sort(max(tmpvoteState(:,stm:end),[],2),'descend');
            data = tmpvoteState(idx,stm:end);
            maxdata = max(data,[],2)/m;
            es = length(find(maxdata>=0.1))+1;
            ss = 2; 
            if ss<es
                mystates{stc} = ss:es;
            else 
                mystates{stc} = ss:ss+1;
            end
            mystates{stc} = 2;
            
            % Converssion happens left to right!
            BI{stc} = de2bi(G{stc}(idx(mystates{stc})));
            DE{stc} = G{stc}(idx(mystates{stc}));
  
        end
        
        HA = {};
        for iii=1:NC
            for jjj = 1:NC
                ci = BI{iii};
                cj = BI{jjj};
                si = size(BI{iii},2);
                sj = size(BI{jjj},2);
                % equalize:
                ci = [ci , zeros(size(BI{iii},1),700-si)];
                cj = [cj , zeros(size(BI{jjj},1),700-sj)];
%                 if si < sj
%                     ci = [ci , zeros(size(BI{iii},1),sj-si)];
%                 else
%                     cj = [cj , zeros(size(BI{jjj},1),si-sj)];
%                 end
                % do all combinations and find mean:
                tmp = [];
                for kk = 1:size(BI{iii},1)
                    for ll = 1:size(BI{jjj},1)
                        tmp = [tmp, pdist2(ci(kk,:),cj(ll,:),'hamming')];
                    end
                end
                if iii == jjj
                    tmp = tmp(find(~eye(size(BI{jjj},1))));
                end
                HA{iii,jjj} = mean(tmp);
            end
        end
        
        img = cell2mat(HA);
        cm = jet(max(max(img))*10000);
        cm(1,:) = [0,0,0];
        
        figure;imagesc(img);
        colormap(cm);
%         title(sprintf('Structured intracluster Hamming (Q=%dms).',Q));
        title(sprintf('Random intracluster Hamming (Q=%dms).',Q));
        
%     end
end

%% confirm not equal states per stimulated cluster:  
for k=1:7
    PBI{k} = padarray(BI{k},[0,700-length(BI{k})],'post');
end
imagesc(cell2mat(PBI'))
%plot single cluster prominent state PLUS stimulation state:
tmpl = zeros(1,700);
tmpl(run.StimMat_str(:,stc)) = 1;
figure;imagesc([PBI{stc};tmpl]');



Qseq = 2:2:50;
for Qs = 1
    stc=1;
    % calculate stimulation overhead:
    
    
    G(idx(mystates))
    
    for k=1:length(mystates)
        
        subplot(ceil(length(mystates)/tria),tria,k);hold on;
        plot(data(mystates(k),:));
        plot(conv(data(mystates(k),:), gaussFilter,'same'),'r', 'linewidth',2);
    end

end


%%
m = run.nruns ; %No of chains (m)
n = tstop; %No of itterations (n)

Qseq = [30];
Nseq = 1;

Qs = 1;
Ns = 1;

Q = Qseq(Qs) ; % simple window (ms)
Qr = floor(n / Q) ; % length of wspiketrain
n_over_b = 20 ; % No of batches
b = floor(Qr/n_over_b); % Batches' length (b)
if b < 2
    error('Batches'' length must be greater than one sample!');
end

s = (((1:Qr)-1) * Q) + 1 ;
e = ((1:Qr) * Q) ;
ts = [];
te = [];
x_len = [];

% Calculate spike trains ONLY for the above selected cells:
wspktrain = zeros(run.nPC,Qr,m);
for ru = 1:m
    spktrain = zeros(run.nPC,n);
    for c=1:run.nPC % must be row vector!!
        if ~isempty(RUNS{1,stc}{c,ru})
            spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
        end
    end
    for k=1:Qr
        wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
    end
end       

% For each Q now, get probability (per run):
%         figure;hold on;
%         cm = hsv(m);
% Find from all the posible states nPC^2, which are represented in my
% sample:
clusterStatesAll = cell(1,7);
voteState = cell(1,7);
for stc=1:NC
    % Calculate spike trains ONLY for the above selected cells:
    wspktrain = zeros(run.nPC,Qr,m);
    for ru = 1:m
        spktrain = zeros(run.nPC,n);
        for c=1:run.nPC % must be row vector!!
            if ~isempty(RUNS{1,stc}{c,ru})
                spktrain(c,round(RUNS{1,stc}{c,ru}.spikes')) = 1;
            end
        end
        for k=1:Qr
            wspktrain(:,k,ru) = any(spktrain(:,s(k):e(k)),2) ;
        end
    end

    % Break to batches to avoid out of memory:
    chunks = 1:700*100:numel(wspktrain);
    clusterStates = {};
    for chunk = 1:length(chunks)-1
        [clusterStates{chunk}, ~, ~] = unique(reshape(wspktrain(chunks(chunk):chunks(chunk+1)-1),run.nPC,[])', 'rows');
    end
    clusterStatesUnion = [];
    clusterStatesUnion = cell2mat(cellfun(@(x)x',clusterStates,'uniformoutput',false));
    [clusterStatesAll{stc}, ~, ~] = unique(clusterStatesUnion', 'rows');
%         % If no batches are used:
%         [clusterStates, randomIA, ~] = unique(reshape(wspktrain,run.nPC,[])', 'rows');
    % indices to each of the unique states in my sample (all runs):
    voteState{1,stc} = zeros(m, size(clusterStatesAll{stc},1));
    for ru = 1:m-1
        ru
            [~,locStates] = ismember(wspktrain(:,:,ru)', clusterStatesAll{stc}, 'rows');
            % Gamhse...
            voteState{1,stc}(ru,:) = accumarray([1:numel(voteState{1,stc}(ru,:)) locStates'].',[voteState{1,stc}(ru,:) ones(1,length(locStates))])' ;
    end
    [~,sidx] = sort(mean(voteState{1,stc}),'descend');
    figure;imagesc(voteState{1,stc}(:,sidx(2:100)))
    figure;plot(sum(voteState{1,stc}(:,sidx(2:100)))/size(clusterStatesAll{stc},1))
end
