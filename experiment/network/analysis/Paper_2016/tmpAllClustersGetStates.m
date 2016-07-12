
mr = -200:5:200; % preference range
% tmpidx = cell(1,length(mr));
% dpsim = cell(1,length(mr));
[X,Y] = meshgrid(1:nPC);
similarity_str = zeros((nPC^2-nPC)/2,3);
similarity_str(:,1) = X(logical(triu(ones(nPC),1)));
similarity_str(:,2) = Y(logical(triu(ones(nPC),1)));
CN_str = m_commonNeighbors(PC2PC_str);
similarity_str(:,3) = CN_str(logical(triu(ones(nPC),1)));
for k=29:length(mr)
    mr(k)
[tmpidx{k},~,dpsim{k},~]=apclusterSparse(similarity_str,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
% [targ_str,~,labels_str] = unique(tmpidx{k}(:,end));

% Random network:
CN_rnd = m_commonNeighbors(PC2PC_rnd);
% tmpidx_rnd = cell(1,length(mr));
% dpsim_rnd = cell(1,length(mr));
similarity_rnd = zeros((nPC^2-nPC)/2,3);
similarity_rnd(:,1) = X(logical(triu(ones(nPC),1)));
similarity_rnd(:,2) = Y(logical(triu(ones(nPC),1)));
similarity_rnd(:,3) = CN_rnd(logical(triu(ones(nPC),1)));
for k=29:length(mr)
    mr(k)
[tmpidx_rnd{k},~,dpsim_rnd{k},~]=apclusterSparse(similarity_rnd,mr(k),'dampfact',0.95, ...
        'convits',200,'maxits',5000,'details');
end
% [targ_rnd,~,labels_rnd] = unique(tmpidx_rnd{k}(:,end));

% figure;plot(cellfun(@(x) x(end), dpsim_rnd))
figure;hold on;
plot(mr,cellfun(@(x) length(unique(x(:,end))),tmpidx_rnd ))
plot(mr,cellfun(@(x) length(unique(x(:,end))),tmpidx ))
legend({'RND','STR'})

%% Array of windows Q to apply:
% Qseq = [100,50,40,30,20,10,8,6,4];
Qseq = [4:2:100];
% genAllClusterStates(run,stc_str,Qseq,tmpidx,batch_str_spikes);
AllClustersStates_str = anlAllClusterStates(run,stc_str,Qseq,700,700,'');

% genAllClusterStates(run,stc_rnd,Qseq,tmpidx_rnd,batch_rnd_spikes);
AllClustersStates_rnd = anlAllClusterStates(run,stc_rnd,Qseq,700,700,'R');
 
% % get unique number of clusters (discard clusters of the same number):
% [~,unc_rnd,~] = unique(cellfun(@(x) length(unique(x(:,end))),tmpidx_rnd));
% 
% for unciters_rnd = 21:length(unc_rnd)
%     [~,~,ulabels_rnd] = unique(tmpidx_rnd{unc_rnd(unciters_rnd)}(:,end))
%     ulabi_rnd = unique(ulabels_rnd);
%     for clustiters_rnd = 302:length(ulabi_rnd)
%         cellsincluster_rnd = (ulabels_rnd == ulabi_rnd(clustiters_rnd));
%         if sum(cellsincluster_rnd) > 1
%             st = batch_str_spikes(cellsincluster_rnd,:);
%             stc=stc_str;
%             configuration = sprintf('NC%dC%d',unciters_rnd,clustiters_rnd);
%             [~, ~] = createVoteState(run, Qseq, st, stc-1, configuration, 'save');
%         end
%     end
% end



