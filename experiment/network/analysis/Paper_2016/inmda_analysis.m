function inmda_analysis(run,inmdaData,inmdaRang,stc,sc,Qseq,configuration)
%% get iNMDA per pyramidal:
% load('C:\\Users\\stefanos\\Documents\\inmdaDetailed_S_PID25.mat');
% inmdaRang = 0:0.0005:0.05;
% inmdaData = inmdaDetailed_S_PID25;
% sc = sc_str;
nmdamean = nan(length(sc),length(sc),50);
wmean = nan(length(sc),length(sc),50);
vmean = nan(length(sc),length(sc),50);
% figure(1);
activation=cell(length(sc),1);
activationgrp=cell(length(sc),1);
activmedian = zeros(length(sc),1);
for jj=1:length(sc)
    totalactiv = [];
    for ii=1:length(sc)
        if ~isempty(inmdaData{sc(ii),sc(jj),1})
            tmp = squeeze(cell2mat( cellfun(@(x) x',inmdaData(sc(ii),sc(jj),:),'uniformoutput',false) ));
            try 
                tmp2 = (inmdaRang(1:end-1)*tmp)./sum(tmp);
            catch e
                disp('Exception: Possibly wrong dims in inmdaRang!')
                tmp2 = [];
            end
            totalactiv = [totalactiv,  (inmdaRang(1:end-1)*tmp)./sum(tmp)];
        end
    end
    activation{jj,1} = totalactiv';
    activationgrp{jj,1} = ones(length(totalactiv),1)*jj;
    activmedian(jj,1) = nanmedian(totalactiv);
end

% x = cell2mat(activation);
% g = cell2mat(activationgrp);
[~,idx] = sort(activmedian,'descend');
newg = [];newx=[];
for jj = 1:length(sc)
    newx = [newx; activation{idx(jj),1}];
    newg = [newg; ones(length(activationgrp{idx(jj),1}),1)*(jj)];
end
figure;boxplot(newx,newg);
title('Synaptic iNMDA per pyramidal.');
xlabel('Pyramidal ID');ylabel('Mean synaptic iNMDA');



%% Get active cells per state:
% Qseq = [10];
% close all;
stateCells = cell(length(Qseq),length(tmp));
localidxStateCells = cell(length(Qseq),length(tmp));
nCoactivePC = cell(length(Qseq),1);
nProminentStatesCheck = 100;
uptocoactive = 50;
for Qi = 1:length(Qseq)
    stateVars = load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab','PID25_iNMDA_Qanalysis_stimulatedClusterOnly_GABAb01_SN2_dvolt',...
        sprintf('cluster_smooth_states_%s_stc%d_SN%d_Q%d_v73.mat',configuration,stc-1,run.sn,Qseq(Qi))));
    delayRange = ceil(1500/Qseq(Qi)):run.tstop/Qseq(Qi) ;
    prominentStates = mean(stateVars.voteState(:,delayRange),2);
    [~,maxfreqidx] = sort(prominentStates,'descend') ;
    nCoactivePC{Qi} = cellfun(@(x) length(regexp(x,'1','match')), stateVars.U(maxfreqidx(1:nProminentStatesCheck),:) );

    tmp = maxfreqidx(1:nProminentStatesCheck);
    for k = 1:length(tmp)
        [~, S] = regexp(stateVars.U{tmp(k)},'1','match');
        stateCells{Qi,k} = sc(S)'; % stimulated cells
        localidxStateCells{Qi,k} = S'; % stimulated cells
    end
end

for Qi = 1:length(Qseq)
    statesPerCell = cell(1,length(sc));

    for k=1:length(sc)
        statesPerCell{k} = find(cell2mat(stateCells(Qi,:)) == sc(k));
    end

    nanidx = isnan(activmedian);
    x=activmedian(~nanidx);
    y=cellfun(@length,statesPerCell)';
    y(nanidx) = [];
    % seperate stimulated that create states from stimulated that don't:
    yzeroidx = find(y==0);
    % yzero = y(yzeroidx);
    xzero = x(yzeroidx);
    y(yzeroidx) = [];
    x(yzeroidx) = [];

    bls = regress(y,[ones(length(x),1) x]);
    [brob,stats] = robustfit(x,y);
    figure;
    scatter(x,y,'filled'); grid on; hold on;
    plot(x,bls(1)+bls(2)*x,'r','LineWidth',2);
    plot(x,brob(1)+brob(2)*x,'g','LineWidth',2)
    legend(sprintf('Cells (%d)',length(x)),'Ordinary Least Squares','Robust Regression')
    title(sprintf('(%s) iNMDA VS state freq (Q=%d)',configuration,Qseq(Qi)));
    xlabel('Median iNMDA (per cell)');ylabel('Freq (active in # states)');

    nonStateRange = linspace(0,max(x),10);
    nonStateHisto = histcounts(xzero,nonStateRange);
    figure;bar(nonStateRange(1:end-1),nonStateHisto)
    title(sprintf('(%s) iNMDA of non state neurons (Q=%d)',configuration,Qseq(Qi)));
    xlabel('Median iNMDA (per cell)');ylabel('Freq of nonactive in states');

    
    coactive_iNMDA = cell(uptocoactive,1);
    groupIdx = cell(uptocoactive,1);
    coactive_no = NaN(uptocoactive,1);
    for k = 1:uptocoactive
        tmpidx = find(nCoactivePC{Qi} == k);
        % means for each synapse for each neuron that participate in the
        % k-th coactive state.
        coactive_iNMDA{k} = cell2mat(activation(unique(cell2mat(localidxStateCells(Qi,tmpidx)')),1));
        groupIdx{k} = ones(length(coactive_iNMDA{k}),1)*k;
        if ~isempty(tmpidx)
            coactive_no(k) = k;
        end
    end
    figure;hold on;
    g = cell2mat(groupIdx);
    x = cell2mat(coactive_iNMDA);
    boxplot(x,g);hold on;
%     scatter(g,x,10,'k','filled');
    coactive_no(isnan(coactive_no)) = [];
    set(gca, 'XTickLabel' ,coactive_no);
    title(sprintf('Q=%d (SN%d)',Qseq(Qi),run.sn));
    xlabel('Co-active PCs');ylabel('iNMDA');
end
% RHO = corr([x,y])
% b1 = x\y;
% yCalc1 = b1*x;
% figure;scatter(x, y);
% hold on;
% plot(x,yCalc1);
