%% Load raw batches:
close all;clear all;
SN = 4;
PID = 25;
for stc=6
%     batch = cell(933,100);
    for ru = 1:50
        pathto = sprintf('/home/cluster/stefanos/Documents/Glia/SAVENMDAupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_r%d_PID%d',stc-1,SN,ru-1,PID);
        batch = load_raw_batch(pathto);
%         inmdaBatch = load_raw_batch(pathto,'inmda');
%         save(sprintf('/home/cluster/stefanos/Documents/Glia/SAVENMDAupdatedStimGABAb01NEWBGST_Ss10c%d_SN%d_inmda_r%d_PID%d.mat',stc-1,SN,ru-1,PID),'batch','-v7.3');
        save(sprintf('/home/cluster/stefanos/Documents/Glia/SAVENMDAupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_mV_r%d_PID%d.mat',stc-1,SN,ru-1,PID),'batch','-v7.3');
        if ~isempty(batch)
            batch_rnd(:,ru) = batch(1:933,:);
        end
    end
end

[~,batch_rnd_spikes] = cellfun(@(x) advanced_spike_count(x,-10,0), batch_rnd, 'uniformoutput', false);
save(sprintf('/home/cluster/stefanos/Documents/Glia/dataParsed2Matlab/SAVENMDAupdatedStimGABAb01NEWBGST_Rs10c%d_SN%d_PID%d_spikes.mat',stc-1,SN,PID),'batch_rnd_spikes','-v7.3');


for ru=1:100
%     ru
%     batch_rnd = batch(:,ru);
        save(sprintf('/home/cluster/stefanos/Documents/Glia/dataParsed2Matlab/updatedStimGABAb01NEWBGST_Rs10c5_SN4_r%d_mV.mat',ru-1),'batch_rnd','-v7.3');

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
configuration = 'rnd';
ru = 20;

eval( sprintf('nc = run.nClusters_%s;',configuration) );
eval( sprintf('cl = clusterLabels_%s;',configuration) );
eval( sprintf('st = batch_%s_spikes;',configuration) );
clusterOrdered = [];
for k=1:nc
    clusterOrdered = [clusterOrdered find(cl==k)'];
end
figure;hold on;
for k=1:700
%     [~,spikes] = advanced_spike_count(batch_str{k,2},-10,0);
    spikes = st{clusterOrdered(k),ru};
    scatter(spikes, ones(1,length(spikes))*k,'.');
end