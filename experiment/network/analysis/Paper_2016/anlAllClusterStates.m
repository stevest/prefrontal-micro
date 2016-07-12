function AllClustersStates = anlAllClusterStates(run,stc,Qseq,maxNC,maxC,erfpref)

% calculate statistics for each cluster/cluster number

% load Qanalysis file:

% get agreed voteStates metric (mean?):

% Compare states across configurations:
close all;
nProminentStatesCheck = 10;
AllClustersStates = cell(maxNC,maxC,length(Qseq));

SN = run.sn;
tstop = run.tstop;
for NC = 1:maxNC
    for C = 1:maxC
        tic;
%         filename = fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('DefaultQanalysis'),...
%                 sprintf('cluster_smooth_states_%sNC%dC%d_stc%d_SN%d_Q%d_v73.mat',erfpref,NC,C,stc-1,SN,Qseq(1)));
                filename = fullfile('C:\Users\stefanos','Documents','DefaultQanalysis',...
                    sprintf('cluster_smooth_states_%sNC%dC%d_stc%d_SN%d_Q%d_v73.mat',erfpref,NC,C,stc-1,SN,Qseq(1)));

        if exist(filename,'file')
            fprintf('NC:%d C:%d ',NC,C);
            parAllClustersStates = cell(1,length(Qseq));
            parfor Qi = 1:length(Qseq)
%                 filename = fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('DefaultQanalysis'),...
%                     sprintf('cluster_smooth_states_%sNC%dC%d_stc%d_SN%d_Q%d_v73.mat',erfpref,NC,C,stc-1,SN,Qseq(Qi)));
                filename = fullfile('C:\Users\stefanos','Documents','DefaultQanalysis',...
                    sprintf('cluster_smooth_states_%sNC%dC%d_stc%d_SN%d_Q%d_v73.mat',erfpref,NC,C,stc-1,SN,Qseq(Qi)));
    %             pairs = AllClustersStates{NC,C,Qi};
                if exist(filename,'file')% && isempty(pairs)
                    loadstruct = load(filename);
                    if  ~isempty(loadstruct.U) && ~isempty(loadstruct.voteState)
                        delayRange = ceil(1500/Qseq(Qi)):tstop/Qseq(Qi) ;
                        meanProminentStates = mean(loadstruct.voteState(:,delayRange),2);
                        medianProminentStates = median(loadstruct.voteState(:,delayRange),2);
                        [~,meanMaxfreqidx] = sort(meanProminentStates,'descend') ;
                        [~,medianMaxfreqidx] = sort(medianProminentStates,'descend') ;
                        nActivePC_str = cellfun(@(x) length(regexp(x,'1','match')), loadstruct.U );
                        nFrequentStates = nProminentStatesCheck;
                        if nFrequentStates > length(loadstruct.U)
                            nFrequentStates = length(loadstruct.U);
                        end
                        parAllClustersStates{Qi} = [nActivePC_str(meanMaxfreqidx(1:nFrequentStates)),...
                            nActivePC_str(medianMaxfreqidx(1:nFrequentStates)),...
                            meanProminentStates(meanMaxfreqidx(1:nFrequentStates)),...
                            medianProminentStates(medianMaxfreqidx(1:nFrequentStates))];
                    end
                end % if filename exist
            end % parfor
            fprintf('parfor %f secs.\n ',toc);
            AllClustersStates(NC,C,:) = parAllClustersStates;
        end
    end
%     tic;
%     save(sprintf('%sNC%d_C%d_AllClustersStates.mat',erfpref,NC,C),'AllClustersStates','-v7.3');
%     fprintf('Saving %f secs.\n',toc);
end
