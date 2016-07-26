function stats = getConfigurationStats(run,p)
eval( sprintf('p.activeCluster = p.stc_%s;',p.activeconfig) );
eval( sprintf('p.stimulatedCells = run.stimulatedCells_%s{p.activeCluster};',p.activeconfig) );
    
analysis = load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',p.postanalysisdir,...
    sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',p.activeconfig,p.activeCluster-1,run.sn,p.Qseq(p.Qi))));
delayRange = ceil(1500/p.Qseq(p.Qi)):run.tstop/p.Qseq(p.Qi);
prominentStates = mean(analysis.voteState(:,delayRange),2);
[maxfreqstates,maxfreqidx] = sort(prominentStates,'descend');
nCoactivePC = cellfun(@(x) length(regexp(x,'1','match')), analysis.U(maxfreqidx(1:p.nProminentStates),:) );

frequencyPerActive = {}; gPC={};
for k=1:p.uptocoactive
    tmp = find(nCoactivePC == k);
    gPC{k,1} = ones(length(tmp),1)*k;
    frequencyPerActive{k,1} = maxfreqstates(tmp);
end

% state intra weight:
stateMeanW = nan(p.nProminentStates,1);
% eval( sprintf('stimulatedCells = sc_%s;',configuration) );
for k=2:p.nProminentStates
    [~, S] = regexp(analysis.U{maxfreqidx(k)},'1','match');
%     length(S) == nCoactivePC(k) ; % SHOULD BE TRUE!
    stateCells = p.stimulatedCells(S)';
    if length(stateCells) > 1
        eval( sprintf('stateCellsW = run.weights_%s(stateCells,stateCells);',p.activeconfig) );
        eval( sprintf('stateCellsC = run.configuration_%s(stateCells,stateCells);',p.activeconfig) );
        % get weights (mean) only for connected pairs:
        stateMeanW(k,1) = mean(stateCellsW(~eye(length(stateCells)) & stateCellsC));
    end
end

% Shuffled state intra weights:
shuffledStateMeanW = nan(p.uptocoactive,1000);
for k=2:p.uptocoactive
    for kk = 1:1000
        shuffledCells = ceil(rand(k,1)*run.nPC)';
        eval( sprintf('shuffledCellsW = run.weights_%s(shuffledCells,shuffledCells);',p.activeconfig) );
        eval( sprintf('shuffledCellsC = run.configuration_%s(shuffledCells,shuffledCells);',p.activeconfig) );
        % get weights (mean) only for connected pairs:
        shuffledStateMeanW(k,kk) = mean(shuffledCellsW(~eye(length(shuffledCells)) & shuffledCellsC));
    end
end

stats.statesFreqTrace = analysis.voteState(maxfreqidx(1:p.nProminentStates),delayRange);
% stats.delayRange = delayRange;
stats.nCoactivePC = nCoactivePC;
% stats.prominentStates = prominentStates;
stats.frequencyPerActive = frequencyPerActive;
stats.gPC = gPC;
% Freuency of only prominent (ie 100) states:
% stats.prominentStates = prominentStates(1:p.nProminentStates);
stats.maxfreqidx = maxfreqidx(1:p.nProminentStates);
stats.maxfreqstates = maxfreqstates(1:p.nProminentStates);
stats.stateMeanW = stateMeanW;
stats.shuffledStateMeanW = shuffledStateMeanW;

end