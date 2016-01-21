clear all;close all;clc;
pathprefix = 'Z:/data/GliaBackup/'
cd(pathprefix)
load('states_07_G.mat')
ST = 7
load(sprintf('%sexperiment_%d\\EXP_ID%d_SN%d_ST%d.mat',pathprefix,12,12,16,ST)) % 8,10

% % Replace connections of PV2PV PC2PV PV2PC with proper ones:
load('C:\Users\marouli\Desktop\matlab_files\PVconnMat_1.mat');
run.state_str(76:100,:,1) = AllConnMat(76:100,:);
run.state_str(:,76:100,1) = AllConnMat(:,76:100);
load('C:\Users\marouli\Desktop\matlab_files\PVconnMat_gapmat.mat');
run.state_str(76:100,76:100,1) = ConnMatPV2PV;
for k = 1: run.nPC
    run.state_str(k,k,1) = 1;
end

% % Replace connections of PV2PV PC2PV PV2PC with proper ones:
load('C:\Users\marouli\Desktop\matlab_files\PVconnMat_1.mat');
run.state_rnd(76:100,:,1) = AllConnMat(76:100,:);
run.state_rnd(:,76:100,1) = AllConnMat(:,76:100);
% load('C:\Users\marouli\Desktop\matlab_files\PVconnMat_gapmat.mat');
run.state_rnd(76:100,76:100,1) = ConnMatPV2PV;
for k = 1: run.nPC
    run.state_rnd(k,k,1) = 1;
end

% run.weights_str(run.weights_str > 1) = 2.34;
% run.exportNetworkParameters(pathprefix)

%%
exportpath = pathprefix;
% Export Network Connectivity:
% Export STRUCTURED parameter matrices in .hoc file:
fprintf('Exporting importGapParametersSTR.hoc...\n');
fid = fopen([exportpath,run.SLASH,'importGapParametersSTR.hoc'],'W');
fprintf(fid,'// This HOC file was generated with MATLAB\n\n');

fprintf(fid,'// Object decleration:\n');
fprintf(fid,'objref gapMatrix\n');
fprintf(fid,'gapMatrix = new Matrix(nPVcells, nPVcells)\n');

fprintf(fid,'\n\n// Import parameters: (long-long text following!)\n\n');
% network connectivity:
for i=1:length(gapmat)
    for j=1:length(gapmat)
        fprintf(fid,'gapMatrix.x[%d][%d] = %d\n',i-1,j-1, gapmat(i,j));
    end
end
fprintf(fid,'//EOF\n');
fclose(fid);

% Export RANDOM parameter matrices in .hoc file:

%%
% run.nPC = nPC;
% run.nPV = nPV;
% close all;
% Each experiment has its own folder now!
run.path = '20secondsRuns_experiment_12_petask' ;
rurange = 1:100;
stop = 20000;
currentCluster = 1;
% load from cluster:
pathprefix = '\\139.91.162.90\cluster\stefanos\Documents\Glia\';
% pathprefix = 'Z:\data\GliaBackup\';
% spathprefix = 'C:\Users\marouli\Documents\';

run.nruns=100;
% RUNS_str = cell(1,run.NC_str);
Sid=1;
fprintf('Loading runs...');
PCcells_str = cell(run.nPC,run.nruns);
% sPCcells_str = cell(run.nPC,run.nruns);
PVcells_str = cell(run.nPC,run.nruns);
% delayPV2PC_str = cell(run.nPC,run.nruns);
% PVcellsEvents = cell(run.nPV,run.nruns);
% PCcellsGABAa_str = cell(run.nPC,run.nruns);
% PCcellsGABAb_str = cell(run.nPC,run.nruns);
% PCcellsica_str = cell(run.nPC,run.nruns);
% PCcellscai_str = cell(run.nPC,run.nruns);
% PCcellsinmda_str = cell(run.nPC,run.nruns);
% PCcellsanyvar_str = cell(run.nPC,run.nruns);
ridx = zeros(run.NC_str(Sid),run.nPC,run.nruns);
basalSegments = 5;
dendseg = {};
failedToLoad = zeros(run.nPC+run.nPV,rurange(end),currentCluster(end));
for stc=currentCluster%:run.NC_str(Sid)
    RUNS_str{1,stc} = cell(run.nPC,run.nruns);
    sRUNS_str{1,stc} = cell(run.nPC,run.nruns);
    for ru = rurange
        fprintf('Run is: %d, of cluster %d\n',ru,stc);
        tic;
        S = cell(1,run.nPC);
        for pc=1:run.nPC
            fprintf('%d,',pc);
            if (run.ISBINARY)
                try
                    %     PCcells_str{pc,ru} = ncell(load(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.txt',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1)),10);
                    if( exist(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1),'file') )
                        PCcells_str{pc,ru} = ncell(nrn_vread(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1),'n'),10);
                        PCcells_str{pc,ru}.clusterID = run.labels_str(pc,Sid);
                        %             sPCcells_str{pc,ru}.clusterID = run.labels_str(pc,Sid);
                        %             [S{pc},~,~] = findUPstates(PCcells_str{pc,ru}.mv(run.stimend*run.dt:run.dt:end),4, 10, -66, 3000 );
                        % Use spikecount if ncell.spike_count is broken:
                        [~,PCcells_str{pc,ru}.spikes] = spike_count(PCcells_str{pc,ru}.mv);
                        PCcells_str{pc,ru}.spikes = single(PCcells_str{pc,ru}.spikes);
                        %             sPCcells_str{pc,ru}.spikes = single(sPCcells_str{pc,ru}.spikes);
                        if (sum(PCcells_str{pc,ru}.spikes>run.stimend)/((stop-run.stimend)/1000)) > 7
                            PCcells_str{pc,ru}.persistentActivity = 1;
                        else
                            PCcells_str{pc,ru}.persistentActivity = 0;
                        end
                        PCcells_str{pc,ru}.mv = [];
                    else
                        failedToLoad(pc,ru,stc) = 1;
                    end
                    %                     PCcellsGABAa_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/GABAa_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     PCcellsGABAb_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/GABAb_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     PCcellsica_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/ica_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     PCcellscai_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/cai_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     PCcellsinmda_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/inmda_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     PCcellsanyvar_str{pc,ru} = load(sprintf('%s%s/RND_SN%d_ST%d/inmda_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
                    %                     sPCcells_str{pc,ru} = ncell(load(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.txt',spathprefix,run.path,run.sn,run.state,stc-1,pc-1,ru-1)),10);
                catch err
                    warning('on','verbose')
                    warning(err.message)
                    ridx(stc,pc,ru) = 1; % true, if file not found
                    continue;
                end
            else
            end
            
            
            
            %             if ~isempty(S{pc});
            %                 PCcells_str{pc,ru}.persistentActivity = 1;
            % %                 sPCcells_str{pc,ru}.persistentActivity = 1;
            %             else
            %                 PCcells_str{pc,ru}.persistentActivity = 0;
            % %                 sPCcells_str{pc,ru}.persistentActivity = 0;
            %             end
            %             PCcells_str{pc,ru} = PCcells_str{pc,ru} ;
            %             sPCcells_str{pc,ru} = sPCcells_str{pc,ru} ;
            % Load delays PV2PC ( for each target)
            %             delayPV2PC_str{pc,ru} = load(sprintf('%s%s/STR_SN%d_ST%d_gaba/delayPV2PC_trg_%d_runs_%d.txt',pathprefix,run.path,run.sn,run.state,pc-1,ru-1));
            
        end
        for pv = nPC+1:nPC+nPV%76:100%76:88
            fprintf('%d,',pv);
            if( exist(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pv-1,ru-1),'file') )
                PVcells_str{pv,ru} = ncell(nrn_vread(sprintf('%s%s/STR_SN%d_ST%d/%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,pv-1,ru-1),'n'),10);
                [~,PVcells_str{pv,ru}.spikes] = spike_count(PVcells_str{pv,ru}.mv);
                PVcells_str{pv,ru}.mv = [];
            else
                failedToLoad(pv,ru,stc) = 1;
            end
        end
%         runtime = toc - tic
    end
    RUNS_str{1,stc} = PCcells_str(:,:);
    %     sRUNS_str{1,stc} = sPCcells_str(:,:);
end
save('\\139.91.162.90\cluster\stefanos\Documents\Glia\RUNS_str_STC0_20S.mat','RUNS_str','-v7.3');
%% Save memory!
clear PCcells_str;
% clear sPCcells_str;

% get vector with valid runs:
find(~any(any(failedToLoad,3)))

fprintf('DONE!\n');

% Plot instantaneous firing frequency:
IFRArray = zeros(nPC,stop);
for k=1:nPC
    spikes = RUNS_str2{1,stc}{k,1}.spikes;
    if ~isempty(spikes)
    IFRArray(k,:) = instantfr(spikes,linspace(0,stop,stop));
    end
end

cm = jet(max(IFRArray(:))*100);
cm(1,:)=[0,0,0];
figure;imagesc(IFRArray(run.StimMat_str(:),:));
% figure;imagesc(IFRArray);
colormap(cm);


% set ONE labels variable:
labels = run.labels_str;
NC = run.NC_str;
StimMat = run.StimMat_str;
cellsPerCluster = run.cellsPerCluster_str;
% Freq and GABA plots

steps = run.dt;
bin = 100;
start = run.stimstart;

bins = start:bin:stop;
freqBinsCell = cell(1,36);
CV = zeros(run.nPC,NC,36);
% Fano = zeros(run.nPC,run.NC_str,36);
ISI = cell(run.nPC,NC,36);
st = start+1:bin:stop-bin+1;
sp = start+bin:bin:stop;
% spikeBins = zeros(run.nPC,5000,36);
spikeBins = struct();
spikeBinsCluster = struct();
spikeBinsNoCluster = struct();
vCell = cell(1,36);
% GABAaCell = cell(1,36);
% GABAbCell = cell(1,36);
% icaCell = cell(1,36);
% caiCell = cell(1,36);
% inmdaCell = cell(1,36);
clstidx = [];
% clstidx = [clstidx ; run.StimMat_str(:,stc)];
% tmp = 1:nPC;
% [ia,ic] = ismember(tmp,run.StimMat_str(:,stc))
% clstidx = [clstidx ; tmp(~ia)'];
% Concatinate stim indices for a clear view of clusters and rest of net:
clstidx = StimMat(:);


% for k = 1:run.NC_rnd
%     clstidx = [clstidx ; find(run.labels_rnd==k)];
% end

if (length(rurange)>1)
    for k=1:run.nPC
        spikeBins(k).spikes = logical(zeros(7,stop));
    end
    clusterCells = find(labels==currentCluster)';
    noClusterCells = find(labels~=currentCluster)';
    for k=1:length(clusterCells)
        spikeBinsCluster(k).spikes = logical(zeros(7,stop));
    end

    for k=1:length(noClusterCells)
        spikeBinsNoCluster(k).spikes = logical(zeros(7,stop));
    end
end

% rurange = [2,4:9]

for ru = rurange
    freqBins = zeros(run.nPC,length(bins),NC);
    inmdaBins = zeros(run.nPC,length(bins),NC);
%     vBins = zeros(run.nPC,length(st),NC);
    for stc = currentCluster%:run.NC_str
        for k = 1:length(clstidx)
            
            spikes = RUNS_str{1,stc}{clstidx(k),ru}.spikes;
            if (length(rurange)>1)
%                 spikeBins(clstidx(k),round(spikes),ru) = 1;
%                 spikes(spikes<start) = [];
                if any(ismember(clusterCells,clstidx(k)))
                    clOvrh = cumsum(cellsPerCluster(1:currentCluster-1));
                    spikeBinsCluster(k-(clOvrh(end))).spikes(ru,round(spikes)) = 1;
                else

                    spikeBinsNoCluster(k).spikes(ru,round(spikes)) = 1;
                end
                spikeBins(clstidx(k)).spikes(ru,round(spikes)) = 1;
            end
            
            freqhisto = (histc(spikes, bins) / (bin/1000));
%             if ~isempty(spikes)
%                 CV(clstidx(k),stc,ru) = std(diff(spikes)) / mean(diff(spikes));
%                 ISI{clstidx(k),stc,ru} = diff(spikes(spikes>1500));
%             end
            if ~isempty(freqhisto)
                freqBins(clstidx(k),:,stc) = freqhisto ;
            end
%             v = RUNS_str{1,stc}{clstidx(k),ru}.mv;
%             v = v(1:steps:end);
%             for kk = 1:length(st)
%                 vBins(clstidx(k),kk,stc) = mean(v(st(kk):sp(kk)));
%             end
%             for kk = 1:length(st)
%                 inmdaBins(clstidx(k),kk,stc) = mean(PCcellsinmda_str{clstidx(k),ru}(st(kk):sp(kk)));
%             end
        end
    end
%     sum(freqBins(:))
    freqBinsCell{ru} = freqBins;
%     inmdaCell{ru} = inmdaBins;
%     vCell{ru} = vBins;
end

% inmda = zeros(nPC,stop);
% GABAa = zeros(nPC,stop);
% GABAb = zeros(nPC,stop);
% ica = zeros(nPC,stop);
% cai = zeros(nPC,stop);
% anyvar = zeros(nPC,stop);
% for k = 1:nPC
%     inmda(k,:) = PCcellsinmda_str{clstidx(k),ru}(1:10:end)';
%     GABAa(k,:) = PCcellsGABAa_str{clstidx(k),ru}(1:10:end)';
%     GABAb(k,:) = PCcellsGABAb_str{clstidx(k),ru}(1:10:end)';
%     ica(k,:) = PCcellsica_str{clstidx(k),ru}(1:10:end)';
%     cai(k,:) = PCcellscai_str{clstidx(k),ru}(1:10:end)';
%     anyvar(k,:) = PCcellscai_str{clstidx(k),ru}(1:10:end)';
% end
% figure('renderer','opengl');
% figure;surf(anyvar(1:40,:),'EdgeColor','none','LineStyle','none');
% figure;surf(inmda(1:40,:),'EdgeColor','none','LineStyle','none');
% figure;surf(Allv(41:100,:));
% figure;surf(Allv(1:40,:));

% %%
% % load PMDdata % PMD dataset from monkey G
% % target onset is at 400 ms
% % all delays > = 400 ms
% 
% fanoParams.alignTime = 0; % this time will become zero time
% fanoParams.boxWidth = 1000; % 50 ms sliding window.
% fanoParams.binSpacing = 2;
% fanoParams.matchReps = 10;
% 
% times = fanoParams.boxWidth:25:stop-fanoParams.boxWidth; 
% % times = 200:25:850; fanoParams.alignTime = 400;fanoParams.boxWidth = 50;% from 200 ms before target onset until 450 ms after.
% Result = VarVsMean(spikeBinsCluster, times, fanoParams);
% plotParams.lengthOfTimeCal = 1000;
% plotParams.leftEdgeOfTimeCal = 500;
% plotParams.plotRawF = 1;
% plotFano(Result,plotParams);
% 
% % scatterParams.axLim = 'auto'; 
% % scatterParams.axLen = 5;
% % plotScatter(Result, -100, scatterParams);
% % text(2.5, 7, '100 ms before target', 'hori', 'center');
% % plotScatter(Result, 100, scatterParams);
% % text(2.5, 7, '100 ms after target', 'hori', 'center');
% % plotScatter(Result, 300, scatterParams);
% % text(2.5, 7, '300 ms after target', 'hori', 'center');
% 
% %%

cl=stc
win = bin;
wst = start+1:win:stop-win+1;
wsp = start+win:win:stop;
cAllFF = cell(5,4);
% cAllv = cell(5,4);
% cAllinmda = cell(5,4);
for ru = rurange
    j = floor((ru-1)/6)+1;
    k = mod(ru-1,6)+1;
    cAllFF{j,k} = freqBinsCell{ru}(clstidx,:,cl);
%     cAllv{j,k} = vCell{ru}(clstidx,:,cl);
%     cAllinmda{j,k} = inmdaCell{ru}(clstidx,:,cl);
end
AllFF = cell2mat(cAllFF);
% Allv = cell2mat(cAllv);
% Allinmda = cell2mat(cAllinmda);

figure();imagesc(AllFF);hold on;
title('Firing Frequency');
set(gca,'YDir','normal');
% for j=1:size(cAllFF,2)
%     plot([j*(length(wst)+1),j*(length(wst)+1)],[0,run.nPC*size(cAllFF,1)],'w','linewidth',2);
%     plot([(j-1)*(length(wst)+1)+run.stimdur/bin,(j-1)*(length(wst)+1)+run.stimdur/bin],[0,run.nPC*size(cAllFF,1)],'r');
% end
% for k=1:size(cAllFF,1)
%     plot([0,size(cAllFF,2)*length(wst)],[k*run.nPC,k*run.nPC],'w');
% end
% set(gca,'Xtick',(length(wst))/2:(length(wst))/2:(length(wst))*size(cAllv,2));
set(gca,'Xticklabel',run.stimend+bin*5 : bin*5 :stop-1);
% set(gca,'Ytick',run.nPC/2:run.nPC:run.nPC*size(cAllFF,1));
% set(gca,'Yticklabel',[1.5,1.9,2.3,2.7,3.0]);
colorbar
%%
% xlabel('GABAb/GABAa');
% ylabel('NMDA/AMPA');

figure();imagesc(Allv);hold on;
title('Voltage')
set(gca,'YDir','normal');
for j=1:size(cAllv,2)
    plot([j*(length(wst)),j*(length(wst))],[0,run.nPC*size(cAllv,1)],'w','linewidth',2);
    plot([(j-1)*(length(wst))+run.stimdur/bin,(j-1)*(length(wst))+run.stimdur/bin],[0,run.nPC*size(cAllv,1)],'r');
end
for k=1:size(cAllv,1)
    plot([0,size(cAllv,2)*length(wst)],[k*run.nPC,k*run.nPC],'w');
end
set(gca,'Xtick',(length(wst))/2:(length(wst)):(length(wst))*size(cAllv,2));
set(gca,'Xticklabel',0.4:0.2:1.4);
set(gca,'Ytick',run.nPC/2:run.nPC:run.nPC*size(cAllv,1));
set(gca,'Yticklabel',[1.5,1.9,2.3,2.7,3.0]);
xlabel('GABAb/GABAa');
ylabel('NMDA/AMPA');


% figure();surf(Allv,'EdgeColor','none','LineStyle','none','FaceLighting','phong');hold on;
% title('Voltage')
% set(gca,'Xtick',(length(wst))/2:(length(wst)):(length(wst))*size(cAllv,2));
% set(gca,'Xticklabel',0.4:0.2:1.4);
% set(gca,'Ytick',run.nPC/2:run.nPC:run.nPC*size(cAllv,1));
% set(gca,'Yticklabel',[1.5,1.9,2.3,2.7,3.0]);
% xlabel('GABAb/GABAa');
% ylabel('NMDA/AMPA');




cells = find(run.labels_str==stc)'
cellslen = length(cells);
figure;hold on;
cm = hsv(cellslen);
tmp2 = [];
tmp = cell(1,cellslen);
histo = cell(1,cellslen);


for k=1:cellslen
    k
    for l = 1:length(ru)
        tmp{k} = [tmp{k}, ISI{cells(k),stc,ru}];
    end
%     tmp2 = [tmp2, tmp];
    histo{k} = histc(tmp{k}, 0:1:60);

%     maxVal(max(histo{k})>maxVal) = max(histo{k});
%     histoconv = conv(histo{k},fspecial('gaussian',[1 10], 2),'same');
%  histo{k} = (histoconv./max(histoconv))*max(histo{k}) ;
     plot(histo{k});
%     plot(conv(histo{k},fspecial('gaussian',[1 10], 2),'same'));
end
% maxVal = max(cell2mat(cellfun(@max, histo, 'Uniformoutput',false)))
% for k=1:cellslen
%     k
%     if (~isempty(tmp{k}))
%     plot(mean(tmp{k}),maxVal,'v','Color',cm(k,1:3));
%     end
% end
% axis([0,50,0,maxVal+10])


% % Plot activity of PC
% figure;
% ru=rurange(end);
% for k=61%22:22
%     k
%     if sum(run.StimMat_str(:,stc)==k)
%         plot(RUNS_str{1,stc}{k,ru}.mv,'r');hold on;
% %         plot(sRUNS_str{1,stc}{k,ru}.mv,'g');hold on;
% %         spikes = RUNS_str{1,stc}{k,ru}.spikes;
% %         scatter(spikes*run.dt,ones(1,length(spikes))*40);
% %         pause();
%     else
%         plot(RUNS_str{1,stc}{k,ru}.mv);hold on;
% %         spikes = RUNS_str{1,stc}{k,ru}.spikes;
% %         scatter(spikes*run.dt,ones(1,length(spikes))*40);
%     end
% %     plot(PCcellsica_str{k,ru}*8000-66,'k');
% %     pause();
% %     cla;
% end
% % legend({'parallel','single'});
% 

% c=11;
% figure;hold on;
% tmp = load(sprintf('%s%s/STR_SN%d_ST%d/anyvar_%d_%d_%d.txt',pathprefix,run.path,run.sn,run.state,stc-1,c-1,ru-1));
% % tmp = 140*nrn_vread(sprintf('%s%s/STR_SN%d_ST%d/anyvar_%d_%d_%d.bin',pathprefix,run.path,run.sn,run.state,stc-1,k-1,ru-1),'n');
% plot(tmp*1000000, 'Linewidth', 3);
% plot(RUNS_str{1,stc}{c,ru}.mv)

% for c = run.StimMat_str(:,stc)'
% %     c=31; % 11    54    34    65    32    67    31
%     figure;
%     plot(RUNS_str{1,stc}{c,ru}.mv,'color',[1,0.0,0.0]);hold on;
%     axesSpikes = gca;
%     axesFF = axes('Position',axesSpikes.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
%     hold(axesFF,'on');
%     ff = 1000./diff(RUNS_str{1,stc}{c,ru}.spikes) ;
%     if ~isempty(ff)
%     plot(axesFF,[0,RUNS_str{1,stc}{c,ru}.spikes(2:end),stop],[0,conv(ff,fspecial('gaussian',[1 10], 3),'same'),0],'color',[0.3,0.3,0.3],'linewidth',3);
%     end
%     title(sprintf('Cell #%d.',c));
% end



figure;hold on;
for c=1:nPC %find(run.labels_str==stc)'
%     if (run.labels_str(c) == stc)
    if (ismember(c,run.StimMat_str(:,stc)))
        plot(RUNS_str{1,stc}{c,ru}.mv,'r');hold on;
    else
        plot(RUNS_str{1,stc}{c,ru}.mv,'b');hold on;
    end
%     c
%     pause;cla;
end
figure;hold on;
for c=nPC+1:nPC+nPV%76:100
plot(PVcells_str{c,ru}.mv(1:10:end));hold on;
% pause;cla;
end
subplot(1,10,1:2);hold on;
tmp = [];
for c=1:nPC
tmp(c) = length(RUNS_str{1,stc}{c,ru}.spikes(RUNS_str{1,stc}{c,ru}.spikes>1500)) / (stop/1000);
end
barh(histc(tmp,1:5:max(tmp)));
set(gca,'XDir','reverse');
axis off;
subplot(1,10,4:10);hold on;
for c=1:nPC
% 1000/mean(diff(PVcells_str{c,ru}.spikes))
plot(c,length(RUNS_str{1,stc}{c,ru}.spikes) / (stop/1000),'*b')
end
figure;hold on;
for c=nPC+1:nPC+nPV%76:100
% 1000/mean(diff(PVcells_str{c,ru}.spikes))
plot(c-nPC,length(PVcells_str{c,ru}.spikes) / (stop/1000),'*r')
end

% plot(PVcells_str{c,ru}.mv(1:10:end),'color',[0.5,0.5,0.5]);hold on;
% axesSpikes = gca;
% axesFF = axes('Position',axesSpikes.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
% hold(axesFF,'on');
% plot(axesFF,[0,PVcells_str{c,ru}.spikes(2:end),stop],[0,conv(1000./diff(PVcells_str{c,ru}.spikes),fspecial('gaussian',[1 10], 3),'same'),0],'color',[0.3,0.3,0.3],'linewidth',3);
% title(sprintf('Cell #%d.',c));

% % Plot activity of PC
% figure;
% ru=rurange(end);
% for k=22%22:22
%     k
%     if sum(run.StimMat_str(:,stc)==k)
%         plot(RUNS_str{1,stc}{k,ru}.mv,'r');hold on;
%         plot(sRUNS_str{1,stc}{k,ru}.mv,'k');hold on;
%         spikes = RUNS_str{1,stc}{k,ru}.spikes;
%         scatter(spikes*run.dt,ones(1,length(spikes))*40);
% %         pause();
%     else
%         plot(RUNS_str{1,stc}{k,ru}.mv);hold on;
%         spikes = RUNS_str{1,stc}{k,ru}.spikes;
%         scatter(spikes*run.dt,ones(1,length(spikes))*40);
%     end
%     plot(PCcellsica_str{k,ru}*8000-66,'k');
% %     pause();
% %     cla;
% end
% legend({'parallel','single'});

%% Parse Connectivity file
close all;clear all;clc;
% parse network connectivity:
% An empty array:
importedarray = [];
% Open the data file to read from it:
fid = fopen( 'C:\Users\marouli\Desktop\importNetworkParametersSTR.hoc', 'r' );
% Check that fid is not -1 (error openning the file).
% read each line of the data in a cell array:
cac = textscan( fid, '%s', 'Delimiter', '\n' );
% size(cac{1},1) must equals the # of rows in your data file.
totalRows = size(cac{1},1);
fprintf('Imported %d rows of data!\n',totalRows);
% Close the file as we don't need it anymore:
fclose( fid );

% for total rows in data:
start = 24;
stop = totalRows-1;

for k=start: stop
%     for k=start:(start+10000)-1
%     fprintf('Parsing data on row %d of %d...\n',k,totalRows);
    currentRow = cac{1}{k,1};
    eachRowElement = strsplit(currentRow, '= ');
    importedarray(k-start+1) = str2num(cell2mat(eachRowElement(2)));
end

mymat = reshape(importedarray,400,400);
figure;imagesc(mymat)

%% data from Konnerth, 2013:
y = [20,0,12,6,4,11,6,5,4,3,5,4,4,5,2,1,3,2,3,0,0,2,0,1,0,0,2,1]';
y = reshape([y,y]',1,[]);
x = 1:length(y);
f = fit(x',y','exp1')
plot(f,x,y)

% create PDF sampler:
expPDF = @(x)f.a*exp(f.b*x)

rang = 0:0.1:length(y);
cumarray = [];
for k=1:100
    tmp = rang(find(histc(ceil(rand(1)*2164),cumsum(expPDF(rang)))));
    if ~isempty(tmp)
        cumarray(k) = tmp;
    end
end

figure;
tmp = hist(cumarray,rang);
plot(tmp);
%%
plot(PCcells_str{22,ru}.mv)

figure;
ru=rurange(end);
for k=1:57
    k
    if sum(run.StimMat_str(:,stc)==k)
        plot(PCcells_str{k,ru}.mv,'r');hold on;
        spikes = PCcells_str{k,ru}.spikes;
        scatter(spikes*10,ones(1,length(spikes))*40);
%         pause();
    else
        plot(PCcells_str{k,ru}.mv);hold on;
        spikes = PCcells_str{k,ru}.spikes;
        scatter(spikes*10,ones(1,length(spikes))*40);
    end
%     plot(PCcellsica_str{k,ru}*8000-66,'k');
    pause();
    cla;
end

figure;
ru=rurange(end);
for k=76:100

        plot(PVcells_str{k,ru}.mv,'r');hold on;
        spikes = PVcells_str{k,ru}.spikes;
        scatter(spikes*10,ones(1,length(spikes))*40);
%         pause();

%     plot(PCcellsica_str{k,ru}*8000-66,'k');
    pause();
    cla;
end

plot(PVcells_str{76,ru}.mv);hold on;

RASTER = zeros(100,4000);
for pc = 1:75
    RASTER(pc,round(PCcells_str{pc,ru}.spikes)) = 1;
end
for pv = 76:100
    RASTER(pv,round(PVcells_str{pv,ru}.spikes)) = 1;
end

imagesc(RASTER)



