%% states_analysis_all_mice.m
% KL PLOTS: Compare states across configurations:
clear all;close all;clc;
Qseq = [100,50,40,30,20,10,8,6,4];
nProminentStatesCheck = 10;
cm = lines(nProminentStatesCheck);
rlowessWidth = 40;
uptocoactive = 10; % plot up to 10 coactive cells.
ConvergeValues_str = zeros(length(Qseq),nProminentStatesCheck,6);
ConvergeValues_rnd = zeros(length(Qseq),nProminentStatesCheck,6);
GoF_str = cell(length(Qseq),nProminentStatesCheck,6);
GoF_rnd = cell(length(Qseq),nProminentStatesCheck,6);
stateW_str = cell(length(Qseq),nProminentStatesCheck,6);
stateW_rnd = cell(length(Qseq),nProminentStatesCheck,6);

SNrange = 2:7;

for mouse = 1:length(SNrange)
    sprintf('mouse %d',SNrange(mouse));
    load(fullfile(osDrive(),'Documents','Glia',sprintf('NetworkCreation_SN%d.mat',SNrange(mouse))));
    %Which cluster is stimulated in each configuration:
    % This is build to load the correct data: BUT NEEDS TO BE CHECKED!!!
    run.tstop = 10000;
    run.nruns = 100;
    if SNrange(mouse) == 2
        stc_rnd = 2;
        stc_str = 4;
        run.tstop = 20000;
        run.nruns = 100;
    elseif SNrange(mouse) == 3
        stc_rnd = 3;
        stc_str = 2;
        run.tstop = 20000;
        run.nruns = 100;
    elseif SNrange(mouse) == 4
        stc_rnd = 6;
        stc_str = 6;
    elseif SNrange(mouse) == 5
        stc_rnd = 5;
        stc_str = 3;
    elseif SNrange(mouse) == 6
        stc_rnd = 5;
        stc_str = 4;
    elseif SNrange(mouse) == 7
        stc_rnd = 6;
        stc_str = 6;
    end
    
    %     load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Rs%dc%d_SN%d_spikes.mat',run.tstop/1000,stc_rnd-1, run.sn)));
    %     load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Ss%dc%d_SN%d_spikes.mat',run.tstop/1000,stc_str-1,run.sn)));
    
    % List of stimulated/non-stimulated cells in each configuration:
    sc_rnd = run.stimulatedCells_rnd{stc_rnd};
    nsc_rnd = find(~ismember(1:700,run.stimulatedCells_rnd{stc_rnd}));
    sc_str = run.stimulatedCells_str{stc_str};
    nsc_str = find(~ismember(1:700,run.stimulatedCells_str{stc_str}));
    
    for Qi = 1:length(Qseq)
        configuration = 'rnd';
        eval( sprintf('stc = stc_%s', configuration) );
        load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('Qanalysis_stimulatedClusterOnly_GABAb01_SN%d',run.sn),...
            sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Qseq(Qi))));
        delayRange = ceil(1500/Qseq(Qi)):run.tstop/Qseq(Qi) ;
        prominentStates = mean(voteState(:,delayRange),2);
        [~,maxfreqidx] = sort(prominentStates,'descend') ;
        nActivePC_rnd = cellfun(@(x) length(regexp(x,'1','match')), U );

        eval( sprintf('stimulatedCells = sc_%s;',configuration) );
        for k=2:nProminentStatesCheck
            %         plot(smooth(voteState(maxfreqidx(k),delayRange)',rlowessWidth,'rlowess'),'color',cm(k,:),'linewidth',2);
            %         plot(voteState(maxfreqidx(k),delayRange)','color',cm(k,:),'linewidth',2);
%             [f,gof] = fit(delayRange',voteState(maxfreqidx(k),delayRange)','exp1')
            %         plot(delayRange,f(delayRange),'color',cm(k,:));
%             ConvergeValues_rnd(Qi,k,mouse) = mean(f(delayRange(end-round(length(delayRange)/10):end)));
%             GoF_rnd{Qi,k,mouse} = gof;
            % state intra weight:
            [~, S] = regexp(U{maxfreqidx(k)},'1','match');
            stateCells = stimulatedCells(S)';
            if length(stateCells) > 1
                eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
                stateW_rnd{Qi,k,mouse} = sort(tmpW(~eye(length(stateCells))),'descend');
            end
        end
     

        configuration = 'str';
        eval( sprintf('stc = stc_%s', configuration) );
        load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('Qanalysis_stimulatedClusterOnly_GABAb01_SN%d',run.sn),...
            sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Qseq(Qi))));
        delayRange = ceil(1500/Qseq(Qi)):run.tstop/Qseq(Qi) ;
        prominentStates = mean(voteState(:,delayRange),2);
        [maxfreqstates,maxfreqidx] = sort(prominentStates,'descend') ;
        nActivePC_str = cellfun(@(x) length(regexp(x,'1','match')), U );

        eval( sprintf('stimulatedCells = sc_%s;',configuration) );
        for k=2:nProminentStatesCheck
            %         plot(smooth(voteState(maxfreqidx(k),delayRange)',rlowessWidth,'rlowess'),'linewidth',2);
            %         plot(voteState(maxfreqidx(k),delayRange)','linewidth',2);
%             [f,gof] = fit(delayRange',voteState(maxfreqidx(k),delayRange)','exp1')
            %         plot(delayRange,f(delayRange),'color',cm(k,:));
%             ConvergeValues_str(Qi,k,mouse) = mean(f(delayRange(end-round(length(delayRange)/10):end)));
%             GoF_str{Qi,k,mouse} = gof;
            % state intra weights:
            [~, S] = regexp(U{maxfreqidx(k)},'1','match');
            stateCells = stimulatedCells(S)';
            if length(stateCells) > 1
                eval( sprintf('tmpW = run.weights_%s(stateCells,stateCells);',configuration) );
                stateW_str{Qi,k,mouse} = sort(tmpW(~eye(length(stateCells))),'descend');
            end
        end
        
    end
    
end % for all SN


weights_range = 0:0.5:7;
w_rnd = cell(length(Qseq),1);
w_str = cell(length(Qseq),1);
for Qi = 1:length(Qseq)
    w_rnd{Qi} = sort(cell2mat(reshape(squeeze(stateW_rnd(Qi,:,:)),[],1)),'ascend');
    w_str{Qi} = sort(cell2mat(reshape(squeeze(stateW_str(Qi,:,:)),[],1)),'ascend');
    figure;hold on;
    title(sprintf('Converge Probabilities Histo (Q=%d)',Qseq(Qi)));
    ylabel('Relative Frequency (all mice pooled)');
    xlabel('Weight (Value)');
    histo_rnd = histcounts(w_rnd{Qi},weights_range);
    histo_str = histcounts(w_str{Qi},weights_range);
    plot(weights_range(2:end-1),histo_rnd(2:end),'b');
    plot(weights_range(2:end-1),histo_str(2:end),'r');
    legend('Random','Structured');
end


convprob_range = 0:0.5:7;
cv_rnd = cell(length(Qseq),1);
cv_str = cell(length(Qseq),1);
for Qi = 1:length(Qseq)
    cv_rnd{Qi} = sort(reshape(squeeze(ConvergeValues_rnd(Qi,:,:))',1,[]),'ascend');
    cv_str{Qi} = sort(reshape(squeeze(ConvergeValues_str(Qi,:,:))',1,[]),'ascend');
    figure;hold on;
    title(sprintf('Converge Probabilities Histo (Q=%d)',Qseq(Qi)));
    ylabel('Relative Frequency (all mice pooled)');
    xlabel('Converge Probability (Value)');
    plot(convprob_range(1:end-1),histcounts(cv_rnd{Qi},convprob_range),'b');
    plot(convprob_range(1:end-1),histcounts(cv_str{Qi},convprob_range),'r');
    legend('Random','Structured');
end

%% Weights of clustered neurons:
Qseq = [100,50,40,30,20,10,8,6,4];
nProminentStatesCheck = 10;
cm = lines(nProminentStatesCheck);
rlowessWidth = 40;
uptocoactive = 10; % plot up to 10 coactive cells.
ConvergeValues_str = zeros(length(Qseq),nProminentStatesCheck,6);
ConvergeValues_rnd = zeros(length(Qseq),nProminentStatesCheck,6);
GoF_str = cell(length(Qseq),nProminentStatesCheck,6);
GoF_rnd = cell(length(Qseq),nProminentStatesCheck,6);
stateW_str = cell(length(Qseq),nProminentStatesCheck,6);
stateW_rnd = cell(length(Qseq),nProminentStatesCheck,6);

SNrange = 2:7;

for mouse = 1:length(SNrange)
    sprintf('mouse %d',SNrange(mouse));
    load(fullfile(osDrive(),'Documents','Glia',sprintf('NetworkCreation_SN%d.mat',SNrange(mouse))));
    %Which cluster is stimulated in each configuration:
    % This is build to load the correct data: BUT NEEDS TO BE CHECKED!!!
    run.tstop = 10000;
    run.nruns = 100;
    if SNrange(mouse) == 2
        stc_rnd = 2;
        stc_str = 4;
        run.tstop = 20000;
        run.nruns = 100;
    elseif SNrange(mouse) == 3
        stc_rnd = 3;
        stc_str = 2;
        run.tstop = 20000;
        run.nruns = 100;
    elseif SNrange(mouse) == 4
        stc_rnd = 6;
        stc_str = 6;
    elseif SNrange(mouse) == 5
        stc_rnd = 5;
        stc_str = 3;
    elseif SNrange(mouse) == 6
        stc_rnd = 5;
        stc_str = 4;
    elseif SNrange(mouse) == 7
        stc_rnd = 6;
        stc_str = 6;
    end

    % List of stimulated/non-stimulated cells in each configuration:
    sc_rnd = run.stimulatedCells_rnd{stc_rnd};
    nsc_rnd = find(~ismember(1:700,run.stimulatedCells_rnd{stc_rnd}));
    sc_str = run.stimulatedCells_str{stc_str};
    nsc_str = find(~ismember(1:700,run.stimulatedCells_str{stc_str}));
    
    cluster_weights_rnd = cell(run.nClusters_rnd,1);
    cluster_weights_str = cell(run.nClusters_str,1);
    for cluster=1:run.nClusters_rnd
        clusterCells = run.stimulatedCells_rnd{cluster};
        cluster_weights_rnd{cluster} = sort(reshape(run.weights_rnd(clusterCells,clusterCells),[],1),'ascend');
    end
    for cluster=1:run.nClusters_str
        clusterCells = run.stimulatedCells_str{cluster};
        cluster_weights_str{cluster} = sort(reshape(run.weights_str(clusterCells,clusterCells),[],1),'ascend');
    end
    
    weights_range = 0:0.1:7;
    figure;hold on;
    title(sprintf('Weights Histo '));
    ylabel('Relative Frequency (all mice pooled)');
    xlabel('Weight (Value)');
    for cluster=1:run.nClusters_rnd
        mean(cluster_weights_rnd{cluster})
        histo_rnd = histcounts(cluster_weights_rnd{cluster},weights_range)./numel(cluster_weights_rnd{cluster});
        plot(weights_range(2:end-1),histo_rnd(2:end),'b');
    end
    for cluster=1:run.nClusters_str
        mean(cluster_weights_str{cluster})
        histo_str = histcounts(cluster_weights_str{cluster},weights_range)./numel(cluster_weights_str{cluster});
        plot(weights_range(2:end-1),histo_str(2:end),'r');
    end
    
end
    
    


