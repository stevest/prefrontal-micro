% calculate statistics for each cluster/cluster number

% load Qanalysis file:

% get agreed voteStates metric (mean?):

% Compare states across configurations:
Qseq = [4:2:100];
VARPID = 25;
close all;
nProminentStatesCheck = 10;
cm = lines(nProminentStatesCheck);
rlowessWidth = 40;
uptocoactive = 10; % plot up to 10 coactive cells.
ConvergeValues_str = zeros(length(Qseq),nProminentStatesCheck);
% ConvergeValues_rnd = zeros(length(Qseq),nProminentStatesCheck);
maxNC = 700;
maxC = 700;
AllClustersStates_str = cell(maxNC,maxC,length(Qseq));
AllClustersStates_rnd = cell(maxNC,maxC,length(Qseq));

for NC = 1:maxNC
    for C = 1:maxC
        for Qi = 1:length(Qseq)
            [NC, C, Qi]
            
            configuration = 'rnd';
            eval( sprintf('stc = stc_%s;', configuration) );
            filename = fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('DefaultQanalysis'),...
                sprintf('cluster_smooth_states_RNC%dC%d_stc%d_SN%d_Q%d_v73.mat',NC,C,stc-1,run.sn,Qseq(Qi)));
            pairs = AllClustersStates_rnd{NC,C,Qi};
            if exist(filename,'file') && isempty(pairs)
                load(filename);
                delayRange = ceil(1500/Qseq(Qi)):run.tstop/Qseq(Qi) ;
                % CHANGED to MEDIAN!!!
                prominentStates = median(voteState(:,delayRange),2);
                [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
                nActivePC_str = cellfun(@(x) length(regexp(x,'1','match')), U );
                
                AllClustersStates_rnd{NC,C,Qi} = [prominentStates(maxfreqidx(1:nProminentStatesCheck)),nActivePC_str(maxfreqidx(1:nProminentStatesCheck))];
                
            end % if filename exist
            
            
            configuration = 'str';
            eval( sprintf('stc = stc_%s;', configuration) );
            filename = fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('DefaultQanalysis'),...
                sprintf('cluster_smooth_states_NC%dC%d_stc%d_SN%d_Q%d_v73.mat',NC,C,stc-1,run.sn,Qseq(Qi)));
            pairs = AllClustersStates{NC,C,Qi};
            if exist(filename,'file') && isempty(pairs)
                load(filename);
                delayRange = ceil(1500/Qseq(Qi)):run.tstop/Qseq(Qi) ;
                % CHANGED to MEDIAN!!!
                prominentStates = median(voteState(:,delayRange),2);
                [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
                nActivePC_str = cellfun(@(x) length(regexp(x,'1','match')), U );
                
                AllClustersStates_str{NC,C,Qi} = [prominentStates(maxfreqidx(1:nProminentStatesCheck)),nActivePC_str(maxfreqidx(1:nProminentStatesCheck))];
                
            end % if filename exist
        end
    end
    
end

%% Plot states frequency as a function of co-active cells:
figure;hold on;
cm = lines(maxNC);
for NC = 1:maxNC
    for C = 1:maxC
        for Qi = 1%:length(Qseq)
            pairs = AllClustersStates{NC,C,Qi};
            if ~isempty(pairs)
                scatter(pairs(:,2),pairs(:,1),NC*90,cm(NC,:));
            end
        end
    end
end


