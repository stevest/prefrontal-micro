%% Load raw batches:
for stc=2
%     batch = cell(933,100);
    for ru = 1
        pathto = sprintf('X:\\Documents\\Glia\\updatedStimGABAb01NEWBGST_Rs20c1_SN2_r%d',ru-1);
        tmp = load_raw_batch(pathto);
        if ~isempty(tmp)
            batch(:,ru) = tmp(1:933,:);
        end
    end
%     save(sprintf('X:\\Documents\\Glia\\dataParsed2Matlab\\updatedStimGABAb01NEWBGST_Rs20c%d_SN3_mV.mat',stc-1),'batch_rnd','-v7.3');
end
%%
% batch_rnd(:,problematicRuns) = cell(933,1);
save(sprintf('X:\\Documents\\Glia\\dataParsed2Matlab\\updatedGABAblowOLDBGST_Rs6c%d_mV.mat',stc-1),'batch_rnd','-v7.3');
%% test cell arrays before deletion
for stc=1:1
    batch_str = cell(933,100);
    load(sprintf('X:\\Documents\\Glia\\dataParsed2Matlab\\GABAblowOLDBGST_Ss6c%d_mV.mat',stc-1));
    
    isemptyarray = cellfun(@isempty,batch_str);
    problematicRuns = find(xor(any(isemptyarray),all(isemptyarray)));
    disp(sprintf('problematicRuns sum is: %d',sum(problematicRuns)));
    if sum(problematicRuns)
        pause();
    end
end
%%
% batch_rnd(:,problematicRuns) = cell(933,1);
save(sprintf('X:\\Documents\\Glia\\dataParsed2Matlab\\updatedGABAblowOLDBGST_Rs6c%d_mV.mat',stc-1),'batch_rnd','-v7.3');
%%
figure;plot(batch_str{152,2})
advanced_spike_count(batch_str{152,2},-10,0,'plot');

%%
clusterOrdered = [];
for k=1:run.nClusters_rnd
    clusterOrdered = [clusterOrdered find(clusterLabels_str==k)'];
end
figure;hold on;
for k=1:700
    [~,spikes] = advanced_spike_count(batch_str{k,2},-10,0);
    scatter(spikes, ones(1,length(spikes))*k,'.');
end