function [voteState, U] = createVoteState(run, Qseq, st, stc, configuration, varargin)

sf = false;
pathto = fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab','DefaultQanalysis');

%check if need to save:
for a=1:length(varargin)
    if strcmp(varargin{a},'save');sf=true;end
end
% if exists a path then get it:
if length(varargin) > 1
    pathto = varargin{2};
end
    

for Q = Qseq
    [voteState, U] = getStates(run, Q, st);
%     smoothed_states = zeros(size(voteState));
%     tic;
%     parfor kk=1:size(voteState,1)
%         smoothed_states(kk,:) = smooth(voteState(kk,:)',100,'moving')'; % SSS changed to moving for speed!
%     end
%     fprintf('Generating smoothed_states took: %fs\n',toc);
    if sf
%         save(sprintf('%s\\cluster_smooth_states_%s_stc%d_SN%d_Q%d_v6.mat',pathto,configuration,stc,run.sn,Q),'voteState','U','smoothed_states','-v6');
        savedfile = fullfile(pathto,sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc,run.sn,Q));
%         save(savedfile,'voteState','U','smoothed_states','-v7.3');
        save(savedfile,'voteState','U','-v7.3');
        fprintf('Saved in: %s\n',savedfile);
    end
end