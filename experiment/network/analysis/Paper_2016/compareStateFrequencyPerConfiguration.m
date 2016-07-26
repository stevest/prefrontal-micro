function compareStateFrequencyPerConfiguration(run,p)
%% Compare states across configurations:
% VARPID = 25;
cm = lines(p.nProminentStates);

fhStateTraces_rnd = cell(length(p.Qseq),1);
fhStateTraces_str = cell(length(p.Qseq),1);
fhStateSorted = cell(length(p.Qseq),1);
fhFreqCoactive_rnd = cell(length(p.Qseq),1);
fhFreqCoactive_str = cell(length(p.Qseq),1);
fhIntraWeights_rnd = cell(length(p.Qseq),1);
fhIntraWeights_str = cell(length(p.Qseq),1);
fhShuffledWeights_rnd = cell(length(p.Qseq),1);
fhShuffledWeights_str = cell(length(p.Qseq),1);
for Qi = 1:length(p.Qseq)
    p.Qi = Qi;
    
    % Stats for Random configuration:
    p.activeconfig = 'rnd';
    eval( sprintf('stats_%s = getConfigurationStats(run,p);',p.activeconfig) );
    
    % plot state traces:
    fhStateTraces_rnd{Qi} = figure;
    hold on;
    % EXCLUDE ZERO STATE:
    for k=2:p.nProminentStates
        plot(smooth(stats_rnd.statesFreqTrace(k,:)',p.rlowessWidth,'rlowess'),'color',cm(k,:),'linewidth',2);
    end
    ylim_rnd = get(gca,'YLim');
    title(sprintf('%s Q=%d (SN%d)',p.activeconfig,p.Qseq(p.Qi),run.sn));
    ylabel('Relative Frequency (%)');xlabel('Time (in Q windows)');
    
    
    % Stats for Structured configuration:
    p.activeconfig = 'str';
    eval( sprintf('stats_%s = getConfigurationStats(run,p);',p.activeconfig) );    
    
    % plot state traces:
    fhStateTraces_str{Qi} = figure;
    hold on;
    % EXCLUDE ZERO STATE:
    for k=2:p.nProminentStates
        plot(smooth(stats_str.statesFreqTrace(k,:)',p.rlowessWidth,'rlowess'),'color',cm(k,:),'linewidth',2);
    end
    ylim_str = get(gca,'YLim');
    title(sprintf('%s Q=%d (SN%d)',p.activeconfig,p.Qseq(p.Qi),run.sn));
    ylabel('Relative Frequency (%)');xlabel('Time (in Q windows)');
    
    %Equalize figure Y axis:
    figure(fhStateTraces_rnd{Qi});set(gca,'YLim',max(ylim_rnd,ylim_str));
    figure(fhStateTraces_str{Qi});set(gca,'YLim',max(ylim_rnd,ylim_str));
    
    % Distribution of 
    fhStateSorted{Qi} = figure;hold on;
    plot(stats_rnd.maxfreqstates,'b');
    plot(stats_str.maxfreqstates,'r');
    title(sprintf('States sorted (SN%d)',run.sn));
    ylabel('State Relative Frequency (%)');xlabel('State id');
    legend('Random','Structured');
    
    %State Frequency per number of coactive:
    fhFreqCoactive_rnd{Qi} = figure;hold on;
    boxplot(cell2mat(stats_rnd.frequencyPerActive),cell2mat(stats_rnd.gPC));
    scatter(cell2mat(stats_rnd.gPC),cell2mat(stats_rnd.frequencyPerActive),10,'k','filled');
    title(sprintf('RND State Frequency per co-active PC (SN%d)',run.sn));
    ylabel('State Relative Frequency (%)');xlabel('Co-active PCs');
    fhFreqCoactive_str{Qi} = figure;hold on;
    boxplot(cell2mat(stats_str.frequencyPerActive),cell2mat(stats_str.gPC));
    scatter(cell2mat(stats_str.gPC),cell2mat(stats_str.frequencyPerActive),10,'k','filled');
    title(sprintf('STR State Frequency per co-active PC (SN%d)',run.sn));
    ylabel('State Relative Frequency (%)');xlabel('Co-active PCs');

    % RANDOM plot intra weights:
    % EXCLUDE STATE WITH ONE ACTIVE PC:
    groupIdx = zeros(p.nProminentStates,1);
    for k=2:p.uptocoactive
        groupIdx(find(stats_rnd.nCoactivePC == k)) = k;
    end
    fhIntraWeights_rnd{Qi} = figure;hold on;
    boxplot(stats_rnd.stateMeanW(groupIdx~=0),groupIdx(groupIdx~=0));
    scatter(groupIdx(groupIdx~=0),stats_rnd.stateMeanW(groupIdx~=0),10,'k','filled');
    ylabel('Weight Value');xlabel('Co-active PCs');
    title(sprintf('RND Weight distribution per Coactive PC (Q=%d,SN%d)',p.Qseq(p.Qi),run.sn));
    ylim_a = get(gca,'YLim');
    
    % plot shuffled intra weights:
    % EXCLUDE STATE WITH ONE ACTIVE PC:
    groupIdx = reshape(repmat(1:p.uptocoactive,1000,[]),[],1);
    fhShuffledWeights_rnd{Qi} = figure;hold on;
    data = reshape(stats_rnd.shuffledStateMeanW,[],1);
    boxplot(data(groupIdx~=1),groupIdx(groupIdx~=1));
    scatter(groupIdx(groupIdx~=1),data(groupIdx~=1),10,'k','filled');
    ylabel('Weight Value');xlabel('Co-active PCs');
    title(sprintf('RND Shuffled Weight distribution per Coactive PC (Q=%d,SN%d)',p.Qseq(p.Qi),run.sn));
    ylim_b = get(gca,'YLim');
    
    %Equalize figure Y axis:
    figure(fhIntraWeights_rnd{Qi});set(gca,'YLim',max(ylim_a,ylim_b));
    figure(fhShuffledWeights_rnd{Qi});set(gca,'YLim',max(ylim_a,ylim_b));
    
    
    % STRUCTURED plot intra weights:
    % EXCLUDE STATE WITH ONE ACTIVE PC:
    groupIdx = zeros(p.nProminentStates,1);
    for k=2:p.uptocoactive
        groupIdx(find(stats_str.nCoactivePC == k)) = k;
    end
    fhIntraWeights_str{Qi} = figure;hold on;
    boxplot(stats_str.stateMeanW(groupIdx~=0),groupIdx(groupIdx~=0));
    scatter(groupIdx(groupIdx~=0),stats_str.stateMeanW(groupIdx~=0),10,'k','filled');
    ylabel('Weight Value');xlabel('Co-active PCs');
    title(sprintf('STR Weight distribution per Coactive PC (Q=%d,SN%d)',p.Qseq(p.Qi),run.sn));
    ylim_a = get(gca,'YLim');
    
    % plot shuffled intra weights:
    % EXCLUDE STATE WITH ONE ACTIVE PC:
    groupIdx = reshape(repmat(1:p.uptocoactive,1000,[]),[],1);
    fhShuffledWeights_str{Qi} = figure;hold on;
    data = reshape(stats_str.shuffledStateMeanW,[],1);
    boxplot(data(groupIdx~=1),groupIdx(groupIdx~=1));
    scatter(groupIdx(groupIdx~=1),data(groupIdx~=1),10,'k','filled');
    ylabel('Weight Value');xlabel('Co-active PCs');
    title(sprintf('STR Shuffled Weight distribution per Coactive PC (Q=%d,SN%d)',p.Qseq(p.Qi),run.sn));
    ylim_b = get(gca,'YLim');
    
    %Equalize figure Y axis:
    figure(fhIntraWeights_str{Qi});set(gca,'YLim',max(ylim_a,ylim_b));
    figure(fhShuffledWeights_str{Qi});set(gca,'YLim',max(ylim_a,ylim_b));
    
end
end
