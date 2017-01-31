function [batch]=load_raw_batch_agnostic(pathto,varargin)
% Load experiment batch..
% ...GIVEN that one run is inside the pathto folder.
% One run is equivalent to one MPI batch job (in its own folder!).

% show all warning messages:
warning('on','all');

r = dir(pathto);
if isempty(r)
    return;
end
rf = {r(~[r.isdir]).name};
% Get only .bin files:
files_v_soma = sort(rf(cellfun(@(x) ~isempty(strfind(x,'v_soma')),rf)));
batch = cell(length(files_v_soma),1);
tic;
for fn = 1:length(files_v_soma)
    filename = fullfile(pathto,files_v_soma{fn});
    vtrace = load(filename);
    % dt equals 0.1:
    vtrace = vtrace(1:10:end);
    batch{fn} = vtrace;
end
fprintf('Time: %f secs\n',toc);
outputfile = fullfile(pathto,'batch.mat');
save(outputfile, 'batch');
% remove files to save space and sanity:
for fn = 1:length(files_v_soma)
    filename = fullfile(pathto,files_v_soma{fn});
    delete(filename);	
end
fprintf('Files v_soma deleted!\n');

end
