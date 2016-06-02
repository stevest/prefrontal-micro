% function to load run results
%blindly execute in a batch order the load commands given.
% to load e.g. a specific current give loadParams. 

function load_batch(runParams, loadParams)

% checks and balances
dirarglen = unique(cellfun(@length,runParams.experimentDirArg));
filearglen = unique(cellfun(@length,runParams.experimentFileArg));

if length(dirarglen) ~=1
    error('experimentDirArg : not consistent arg number!');
end

if length(filearglen) ~=1
    error('experimentFileArg : not consistent arg number!');
end

if length(strfind(runParams.experimentDirStr,'%')) ~= dirarglen
    error('experimentDirStr : dir string wrong arg number!');
end

if length(strfind(runParams.experimentFileStr,'%')) ~= filearglen
    error('experimentDirStr : file string wrong arg number!');
end

if length(strfind(runParams.experimentSpikesStr,'%')) ~= length(strfind(runParams.experimentDirStr,'%'))
    error('experimentSpikesStr : str string wrong arg number!');
end

if dirarglen ~= filearglen
    error('arguments : different arg number!');
end

N = length(runParams.experimentDirArg);

if ~isempty(loadParams)
    specfilearglen = unique(cellfun(@length,loadParams.specificsFileArg));
    if length(specfilearglen) ~=1
        error('specificsFileArg : not consistent arg number!');
    end
    if length(strfind(loadParams.specificsFileStr,'%')) ~= specfilearglen
        error('specificsFileStr : file string wrong arg number!');
    end
    N = length(loadParams.specificsFileArg);
end


for ru = 1:N
    pathto = fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentDirStr, runParams.experimentDirArg{ru}{:}));
    if isempty(loadParams)
        batch = load_raw_batch(pathto);
        save(fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentFileStr, runParams.experimentFileArg{ru}{:})),'batch','-v7.3');
        if ~isempty(batch)
            batch_all(:,ru) = batch(1:933,:);
        end
    else
        specificBatch = load_raw_batch(pathto,loadParams.specifics);
        save(fullfile(osDrive(),'Documents','Glia',sprintf(loadParams.specificsFileStr, loadParams.specificsFileArg{ru}{:})),'specificBatch','-v7.3');
    end
end

if isempty(loadParams)
    spikevar = sprintf('batch_%s_spikes',runParams.config);
    eval( sprintf('[~,%s] = cellfun(@(x) advanced_spike_count(x,-10,0), batch_all, ''uniformoutput'', false);',spikevar ))
    save(fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentSpikesStr, runParams.experimentFileArg{ru}{:})),spikevar,'-v7.3');
end

return;
