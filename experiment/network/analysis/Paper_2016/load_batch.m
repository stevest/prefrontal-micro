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

% because spikes str lacks run var:
% if length(strfind(runParams.experimentSpikesStr,'%')) ~= length(strfind(runParams.experimentDirStr,'%'))
%     error('experimentSpikesStr : str string wrong arg number!');
% end

if dirarglen ~= filearglen
    error('arguments : different arg number!');
end

N = length(runParams.experimentDirArg);
batch_all = cell(1,N);

if ~isempty(loadParams)
    specfilearglen = unique(cellfun(@length,loadParams.specificsFileArg));
    if length(specfilearglen) ~=1
        error('specificsFileArg : not consistent arg number!');
    end
    if length(strfind(loadParams.specificsFileStr,'%')) ~= specfilearglen
        error('specificsFileStr : file string wrong arg number!');
    end
    N = length(loadParams.specificsFileArg);
    specificBatch_all = cell(1,N);
end


% load in parallel
tic;
experimentDirStr = runParams.experimentDirStr;
experimentDirArg = runParams.experimentDirArg;
if isempty(loadParams)
    parfor ru = 1:N
        pathto = fullfile(osDrive(),'Documents','Glia',sprintf(experimentDirStr, experimentDirArg{ru}{:}));
        batch = load_raw_batch(pathto);
%         if ~isempty(batch)
            batch_all(1,ru) = {batch(1:933,:)};
%         else
%             warning('Empty data loaded!!!');
%         end
    end
else
    specifics = loadParams.specifics;
    for ru = 1:N
        pathto = fullfile(osDrive(),'Documents','Glia',sprintf(experimentDirStr, experimentDirArg{ru}{:}));
        specificBatch = load_raw_batch(pathto,specifics);

%         if ~isempty(specificBatch)
            specificBatch_all(1,ru) = {specificBatch};
%         else
%             warning('Empty data loaded!!!');
%         end
    end
end
disp(sprintf('Data loaded in %f seconds.',toc));


% % load in parallel
% tic;
% parfor ru = 1:N
%     pathto = fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentDirStr, runParams.experimentDirArg{ru}{:}));
%     if isempty(loadParams)
%         batch = load_raw_batch(pathto);
%         if ~isempty(batch)
%             batch_all(:,ru) = batch(1:933,:);
%         else
%             warning('Empty data loaded!!!');
%         end
%     else
%         specificBatch = load_raw_batch(pathto,loadParams.specifics);
%         
%         if ~isempty(specificBatch)
%             specificBatch_all(:,ru) = specificBatch;
%         else
%             warning('Empty data loaded!!!');
%         end
%     end
% end
% disp(sprintf('Data loaded in %f seconds.',toc));


% save in serial:
for ru = 1:N
    if isempty(loadParams)
        batch = batch_all{1,ru};
        save(fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentFileStr, runParams.experimentFileArg{ru}{:})),'batch','-v7.3');
    else
        specificBatch = specificBatch_all{1,ru};
        save(fullfile(osDrive(),'Documents','Glia',sprintf(loadParams.specificsFileStr, loadParams.specificsFileArg{ru}{:})),'specificBatch','-v7.3');
    end
end

if isempty(loadParams)
    spikevar = sprintf('batch_%s_spikes',runParams.config);
    eval( sprintf('[~,%s] = cellfun(@(x) advanced_spike_count(x,-10,0), batch_all, ''uniformoutput'', false);',spikevar ))
    save(fullfile(osDrive(),'Documents','Glia',sprintf(runParams.experimentSpikesStr, runParams.experimentSpikesArg{:})),spikevar,'-v7.3');
end

return;
