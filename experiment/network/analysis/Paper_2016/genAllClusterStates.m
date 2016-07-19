function genAllClusterStates(run,stc,Qseq,APallidx,AllSpiketrains)

% get unique number of clusters (discard clusters of the same number):
[~,unc,~] = unique(cellfun(@(x) length(unique(x(:,end))),APallidx));

for unciters = 21:length(unc)
    [~,~,ulabels] = unique(APallidx{unc(unciters)}(:,end));
    ulabi = unique(ulabels);
    for clustiters = 1:length(ulabi)
        cellsincluster = (ulabels == ulabi(clustiters));
        if sum(cellsincluster) > 1
            st = AllSpiketrains(cellsincluster,:);
            configuration = sprintf('RNC%dC%d',unciters,clustiters);
            [~, ~] = createVoteState(run, Qseq, st, stc-1, configuration, 'save');
        end
    end
end