% load inmda data:

pathto = sprintf('X:\\Documents\\Glia\\updatedStimiNMDAt10GABAb02NEWBGST_Ss20c1_SN3_r0');

ru=1
tmp = load_raw_batch([pathto,'\voltage']);
if ~isempty(tmp)
    batch_t10(:,ru) = tmp(1:933,:);
end

r = dir(pathto);
rf = {r(~[r.isdir]).name};
% Get only .bin files:
files = sort(rf(cellfun(@(x) strcmp(x(end-3:end),'.bin'),rf)));
files = sort(files(cellfun(@(x) strcmp(x(1:5),'inmda'),files)));
n = size(files,2);

if (n == 0)
   warning('No files found in folder!');
end

inmda = cell(n,1);
parfor fn=1:n
    filename = fullfile(pathto,files{fn});
%     tmp = strsplit(files{fn}, {'srcPC','trgPC','_','.'});
%     sid = str2double(tmp{2})+1;
%     tid = str2double(tmp{3})+1;
%     disp(filename);
    try
        current = nrn_vread(filename,'n');
    catch e
        warning('Error with nrn_vread() !');
    end
    % dt equals 0.1:
    inmda{fn} = current(1:10:end);
end


%% Disperse current array to each connected pair:
coords = zeros(n,2);
for k=1:n
    tmp = strsplit(files{k}, {'srcPC','trgPC','_','.'});
    sid = str2double(tmp{2})+1;
    tid = str2double(tmp{3})+1;
    coords(k,1:2) = [sid, tid];
end

trgs = unique(coords(:,1));
samecell = {};
for k=1:length(trgs)
    samecell{k} = find(coords(:,1)==trgs(k));
end

%%
figure;hold on;
for k=1:length(trgs)
    plot(batch{trgs(k)});
    for j=1:length(samecell{k})
        plot(inmda{samecell{k}(j)});
    end
    pause;
    cla;
end

%%
figure;
for k=1:n
    k
    plot(inmda_str{k});
    pause;
    cla;
end
%%
clusterOrdered = [];
for k=1:run.nClusters_rnd
    clusterOrdered = [clusterOrdered find(clusterLabels_str==k)'];
end
figure;hold on;
for k=1:700
    [~,spikes] = advanced_spike_count(batch_t10{clusterOrdered(k),1},-10,0);
    scatter(spikes, ones(1,length(spikes))*k,'.');
end