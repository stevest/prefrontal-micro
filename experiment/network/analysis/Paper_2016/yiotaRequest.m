%% pyramidal CV ONLY FOR STIMULATED CLUSTER AND DELAY PERIOD!
load(fullfile(osDrive(),'Documents','Glia','NetworkCreation_SN3.mat'));
%Which cluster is stimulated in each configuration:
stc_rnd = 3;
stc_str = 2;
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Rs20c%d_SN%d_spikes.mat',stc_rnd-1, run.sn)));
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Ss20c%d_SN%d_spikes.mat',stc_str-1,run.sn)));
run.tstop = 20000;
run.nruns = 100;

% List of stimulated/non-stimulated cells in each configuration:
sc_rnd = run.stimulatedCells_rnd{stc_rnd};
nsc_rnd = find(~ismember(1:700,run.stimulatedCells_rnd{stc_rnd}));
sc_str = run.stimulatedCells_str{stc_str};
nsc_str = find(~ismember(1:700,run.stimulatedCells_str{stc_str}));

%% Overview: FF and CV:

FF_stim_rnd = zeros(length(sc_rnd),size(batch_rnd_spikes,2));
FF_delay_rnd = zeros(length(sc_rnd),size(batch_rnd_spikes,2));
CV_stim_rnd = zeros(length(sc_rnd),size(batch_rnd_spikes,2));
CV_delay_rnd = zeros(length(sc_rnd),size(batch_rnd_spikes,2));
FF_stim_str = zeros(length(sc_str),size(batch_str_spikes,2));
FF_delay_str = zeros(length(sc_str),size(batch_str_spikes,2));
CV_stim_str = zeros(length(sc_str),size(batch_str_spikes,2));
CV_delay_str = zeros(length(sc_str),size(batch_str_spikes,2));
for ru=1:size(batch_rnd_spikes,2)
    for c=1:length(sc_rnd)
%             [~,st_rnd] = advanced_spike_count(batch_rnd{sc_rnd(c),ru}(1500:end), -10, 0);
        st_rnd = batch_rnd_spikes{sc_rnd(c),ru};
        st_stim_rnd = st_rnd(st_rnd<=1500);
        st_delay_rnd = st_rnd(st_rnd>1500);
        ISI_stim_rnd = diff(st_stim_rnd);
        ISI_delay_rnd = diff(st_delay_rnd);
        FF_stim_rnd(c,ru) = length(st_stim_rnd)/1.5;
        FF_delay_rnd(c,ru) = length(st_delay_rnd)/((run.tstop-1500)/1000);
        CV_stim_rnd(c,ru) = std(ISI_stim_rnd)/mean(ISI_stim_rnd);
        CV_delay_rnd(c,ru) = std(ISI_delay_rnd)/mean(ISI_delay_rnd);
    end
end
for ru=1:size(batch_str_spikes,2)
    for c=1:length(sc_str)
%         [~,st_str] = advanced_spike_count(batch_str{sc_str(c),ru}(1500:end), -10, 0);
        st_str = batch_str_spikes{sc_str(c),ru};
        st_stim_str = st_str(st_str<=1500);
        st_delay_str = st_str(st_str>1500);
        ISI_stim_str= diff(st_stim_str);
        ISI_delay_str = diff(st_delay_str);
        FF_stim_str(c,ru) = length(st_stim_str)/1.5;
        FF_delay_str(c,ru) = length(st_delay_str)/((run.tstop-1500)/1000);
        CV_stim_str(c,ru) = std(ISI_stim_str)/mean(ISI_stim_str);
        CV_delay_str(c,ru) = std(ISI_delay_str)/mean(ISI_delay_str);
    end
end

% Overal delay period activity (FF average per cell):

cmmax = max(max([FF_delay_rnd(:);FF_delay_str(:)]));
cmmin = min(min([FF_delay_rnd(:);FF_delay_str(:)]));

figure;imagesc(FF_delay_rnd);
caxis manual
caxis([cmmin cmmax]);
cm = parula(1000);
cm(1,:) = [0,0,0];
colormap(cm);title('Random');
xlabel('Runs');ylabel('Cell ID');title('FF average per cell (delay period)');
figure;imagesc(FF_delay_str);
caxis manual
caxis([cmmin cmmax]);
colormap(cm);title('Structured');
xlabel('Runs');ylabel('Cell ID');title('FF average per cell (delay period)');

% CV for stimulated cluster (delay period):
edges = 0:0.1:2;
tmp_delay = [histc(CV_delay_rnd(:),edges)'/sum(~isnan(CV_delay_rnd(:)));histc(CV_delay_str(:),edges)'/sum(~isnan(CV_delay_str(:)))];
figure;hold on;plot(edges(2:end),tmp_delay(:,2:end));
cm = get(groot,'DefaultAxesColorOrder');
tmp_stim = [histc(CV_stim_rnd(:),edges)'/sum(~isnan(CV_stim_rnd(:)));histc(CV_stim_str(:),edges)'/sum(~isnan(CV_stim_str(:)))];
plot(edges(2:end),tmp_stim(1,2:end),'color',cm(1,:));
plot(edges(2:end),tmp_stim(2,2:end),'color',cm(2,:));

dmax=max([tmp_delay(:);tmp_stim(:)]);dmax=dmax+dmax*0.1;
scatter(nanmean([CV_delay_rnd(:),CV_delay_str(:)]),[1,1]*dmax,[],cm(1:2,:),'v','fill');
scatter(nanmean([CV_stim_rnd(:),CV_stim_str(:)]),[1,1]*dmax,[],cm(1:2,:),'v','fill');

legend({'Random','Structured'});xlabel('Coefficient of Variation (CV)');ylabel('Count (#)');title('Histogram of CV (Cluster in delay period)');

% FF for stimulated cluster (delay period):
edges = 0:2:60;
tmp = [histc(FF_delay_rnd(:),edges)'/sum(~isnan(FF_delay_rnd(:)));histc(FF_delay_str(:),edges)'/sum(~isnan(FF_delay_str(:)))];
figure;plot(edges(2:end),tmp(:,2:end));
legend({'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Count (#)');title('Histogram of FF (Cluster in delay period)');
set(gca,'XScale','log');


%% Recruitement:
FFr_rnd = zeros(length(nsc_rnd),run.nruns);
FFr_str = zeros(length(nsc_str),run.nruns);

for ru=1:size(batch_rnd_spikes,2)
    for c=1:length(nsc_rnd)
        st_rnd = batch_rnd_spikes{nsc_rnd(c),ru};
        st_rnd = st_rnd(st_rnd>1500);
        FFr_rnd(c,ru) = length(st_rnd)/((run.tstop-1500)/1000);
    end
end
for ru=1:size(batch_str_spikes,2)
    for c=1:length(nsc_str)
        st_str = batch_str_spikes{nsc_str(c),ru};
        st_str = st_str(st_str>1500);
        FFr_str(c,ru) = length(st_str)/((run.tstop-1500)/1000);
    end
end

edges = 0:0.2:8;
tmp = [histc(FFr_rnd(:),edges)'/sum(~isnan(FFr_rnd(:)));histc(FFr_str(:),edges)'/sum(~isnan(FFr_str(:)))];
% figure;plot(edges,histcounts(FFr_rnd(FFr_rnd~=0),edges))
% figure;plot(edges,histcounts(FFr_str(FFr_str~=0),edges))
figure;hold on;plot(edges,tmp)
legend({'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Relative Frequency');title('Non-Stimulated Cells Firing Rate');
set(gca,'XScale','log');
plot([0.4,0.4],[0,0.04])

% Sort by mean activity (per cell):
[MEAN_rnd,IDX_rnd] = sort(mean(FFr_rnd,2),'descend');
[MEAN_str,IDX_str] = sort(mean(FFr_str,2),'descend');

% Set limit in background activity (1Hz):
figure;hold on;
plot(MEAN_rnd(1:40));
plot(MEAN_str(1:40),'r');
xlabel('Cell ID (#)');ylabel('Mean Firing Frequency');

cmmax = max(max([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,:)]));
cmmin = min(min([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,:)]));

figure;imagesc(FFr_rnd(IDX_rnd,:));
caxis manual
caxis([cmmin cmmax]);
cm = hot(1000);
cm(1,:) = [0,0,0];
colormap(cm);title('Random');
figure;imagesc(FFr_str(IDX_str,:));
caxis manual
caxis([cmmin cmmax]);
colormap(cm);title('Structured');



%% Load Random
for k=1:50
    pathto = sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\GABAb02NEWBGST_Rs20c0r%d',k-1);
    if exist(pathto,'dir')
        [batch_rnd(:,k),tstop] = load_raw_batch(pathto);
    end
end
%% Load Structured
for k=48:80
    pathto = sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\GABAb02NEWBGST_Ss20c0r%d',k-1);
    if exist(pathto,'dir')
        [batch_str(:,k),tstop] = load_raw_batch(pathto);
    end
end
%% Load Background
batch_bg = {};
for k=1
    pathto = sprintf('\\\\139.91.162.90\\cluster\\stefanos\\Documents\\Glia\\BG_Ss20c0r%d',k-1);
    if exist(pathto,'dir')
        [batch_bg(:,k),tstop] = load_raw_batch(pathto);
    end
end
%% Find Up States ONLY FOR STIMULATED CLUSTER AND DELAY PERIOD!
stc = 2;
sc_rnd = run.stimulatedCells_rnd{stc};
NumMat_rnd = zeros(100,1);
UPFFMat_rnd = cell(100,1);
DurMat_rnd = cell(100,1);

% edges = -80:2:10;
% vth = {};
for ru=1:80
    for c=1:length(sc_rnd)
        [~,DurMat_rnd{c,ru},UPs_rnd]=find_up_states(batch_rnd{sc_rnd(c),ru}(1500:end),2,18,-66,0);
        NumMat_rnd(c,ru) = length(DurMat_rnd{c,ru});
        UPFFMat_rnd{c,ru} = (cellfun(@(x) advanced_spike_count(x, -10, 0),UPs_rnd)/18.5)';
%         vth{c,ru} = histcounts(batch_rnd{sc_rnd(c),ru}(1500:end),edges);
%         vth{c,ru} = mean(batch_rnd{sc_rnd(c),ru}(1500:end));
% vth{c,ru} = mean(batch_str{sc_str(c),ru}(1500:end));
% vth{c,ru} = histcounts(batch_str{sc_str(c),ru}(1500:end),edges);
    end
end
% figure;hold on;
% for c=1:100
% plot(vth{c,2});
% % scatter(vth{c,2},1,'+');
% % pause();cla;
% end
% title('Structured');

sc_str = run.StimMat_str(:,1);
NumMat_str = zeros(100,1);
UPFFMat_str = cell(100,1);
DurMat_str = cell(100,1);

for ru=1:80
    for c=1:length(sc_str)
        [~,DurMat_str{c,ru},UPs_str]=find_up_states(batch_str{sc_str(c),ru}(1500:end),2,18,-66,0);
        NumMat_str(c,ru) = length(DurMat_str{c,ru});
        UPFFMat_str{c,ru} = (cellfun(@(x) advanced_spike_count(x, -10, 0),UPs_str)/18.5)';
    end
end

%FF of spikes uppon Up states
edges = 0:0.2:50;
figure;plot(edges(1:end-1),histcounts(cell2mat(reshape(UPFFMat_rnd,1,[])),edges)/80);
hold on;plot(edges(1:end-1),histcounts(cell2mat(reshape(UPFFMat_str,1,[])),edges)/80,'r');
legend({'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Count (#)');title('Histogram FF uppon UP state');
%Duration of Up states
edges = 0:0.1:20;
figure;plot(edges(1:end-1),histcounts(cell2mat(reshape(DurMat_rnd,1,[]))/1000,edges)/80);
hold on;plot(edges(1:end-1),histcounts(cell2mat(reshape(DurMat_str,1,[]))/1000,edges)/80,'r');
legend({'Random','Structured'});xlabel('Duration in seconds');ylabel('Count (#)');title('Duration of UP state');
% Number of UP states:
edges = 0:1:5;
tmp = [histcounts(NumMat_rnd(:),edges)',histcounts(NumMat_str(:),edges)']/80;
figure;bar(edges(2:end-1),tmp(2:end,:));
legend({'Random','Structured'});xlabel('No of Up states ');ylabel('Count (#)');title('Number of UP state');




%% Compute Synchronicity (a big issue!):
% SPIKY_loop_results_cluster_rnd = cell(1,72);
% SPIKY_loop_results_cluster_str = cell(1,72);
% SPIKY_loop_results_network_rnd = cell(1,72);
% SPIKY_loop_results_network_str = cell(1,72);
SPIKY_loop_results_recruit_rnd = cell(1,72);
SPIKY_loop_results_recruit_str = cell(1,72);
for ru = 2:72
    ru
        % construct input array:
%         spikes_cluster_rnd = batch_rnd_spikes(sc_rnd(:),ru);
%         spikes_cluster_str = batch_str_spikes(sc_str(:),ru);
%         spikes_network_rnd = batch_rnd_spikes(1:700,ru);
%         spikes_network_str = batch_str_spikes(1:700,ru);
        spikes_recruit_rnd = batch_rnd_spikes(nsc_rnd(:),ru);
        spikes_recruit_str = batch_str_spikes(nsc_str(:),ru);
        %Construct parameters struct:
        para.tmin = 0;
        para.tmax = 20000;
        para.dts = 1;
        
        m_para.all_measures_string = {'SPIKE_Synchro';'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_forward';'PSTH';};  % order of select_measures
        para.select_measures = [1 0 0 0 0 0];  % Select measures (0-calculate,1-do not calculate)
        %Compute synchronicity:
%         para.num_trains = length(spikes_cluster_rnd);
%         try
%         SPIKY_loop_results_cluster_rnd{ru} = SPIKY_loop_f_distances(spikes_cluster_rnd,para);
%         catch e
%             SPIKY_loop_results_cluster_rnd{ru}.SPIKE_Synchro = [];
%         end
%         SPIKY_loop_results_cluster_rnd{ru}.Spikes = [];
% %         SPIKE_profile_rnd{1,ru} = SPIKY_loop_results_rnd.SPIKE.profile;
        
%         para.num_trains = length(spikes_cluster_str);
%         try
%         SPIKY_loop_results_cluster_str{ru} = SPIKY_loop_f_distances(spikes_cluster_str,para);
%         catch e
%             
%             SPIKY_loop_results_cluster_str{ru}.SPIKE_Synchro = [];
%         end
%         SPIKY_loop_results_cluster_str{ru}.Spikes = [];
% %         SPIKE_profile_str{1,ru} = SPIKY_loop_results_str.SPIKE.profile;

%         para.num_trains = length(spikes_network_rnd);
%         SPIKY_loop_results_network_rnd{ru} = SPIKY_loop_f_distances(spikes_network_rnd,para);
%         SPIKY_loop_results_network_rnd{ru}.Spikes = [];
%         
%         para.num_trains = length(spikes_network_str);
%         SPIKY_loop_results_network_str{ru} = SPIKY_loop_f_distances(spikes_network_str,para);
%         SPIKY_loop_results_network_str{ru}.Spikes = [];
        
        para.num_trains = length(spikes_recruit_rnd);
        try
        SPIKY_loop_results_recruit_rnd{ru} = SPIKY_loop_f_distances(spikes_recruit_rnd,para);
        catch e
            SPIKY_loop_results_recruit_rnd{ru}.SPIKE_Synchro = [];
            disp('Xtypisa')
        end
        SPIKY_loop_results_recruit_rnd{ru}.Spikes = [];
        
        para.num_trains = length(spikes_recruit_str);
        try
        SPIKY_loop_results_recruit_str{ru} = SPIKY_loop_f_distances(spikes_recruit_str,para);
        catch e
            SPIKY_loop_results_recruit_str{ru}.SPIKE_Synchro = [];
            disp('Xtypisa')
        end
        SPIKY_loop_results_recruit_str{ru}.Spikes = [];
end

% % Get smoothed  SPIKE_synchronous:
% for k = 1:80%length(SPIKY_loop_results_rnd)
%     SPIKY_loop_results_rnd{k}.SPIKE_Synchro.sprofile = smooth(SPIKY_loop_results_rnd{k}.SPIKE_Synchro.profile,100,'rlowess');
% end
% for k = 1:80%length(SPIKY_loop_results_str)
%     SPIKY_loop_results_str{k}.SPIKE_Synchro.sprofile = smooth(SPIKY_loop_results_str{k}.SPIKE_Synchro.profile,100,'rlowess');
% end

irange = 1:20:19999;
SPIKY_interp_rnd = zeros(length(SPIKY_loop_results_rnd),length(irange));
SPIKY_interp_str = zeros(length(SPIKY_loop_results_str),length(irange));


for k = 1:length(SPIKY_loop_results_rnd)
    [~,iidx,~] = unique(SPIKY_loop_results_rnd{k}.SPIKE_Synchro.time);
    SPIKY_interp_rnd(k,:) = interp1(SPIKY_loop_results_rnd{k}.SPIKE_Synchro.time(iidx),SPIKY_loop_results_rnd{k}.SPIKE_Synchro.profile(iidx),irange);
%     hr = plot(SPIKY_loop_results_rnd{k}.SPIKE_Synchro.time,SPIKY_loop_results_rnd{k}.SPIKE_Synchro.sprofile,'color',cm(1,:));
end
    
for k = 1:length(SPIKY_loop_results_str)
    [~,iidx,~] = unique(SPIKY_loop_results_str{k}.SPIKE_Synchro.time);
    SPIKY_interp_str(k,:) = interp1(SPIKY_loop_results_str{k}.SPIKE_Synchro.time(iidx),SPIKY_loop_results_str{k}.SPIKE_Synchro.profile(iidx),irange);
%     hs = plot(SPIKY_loop_results_str{k}.SPIKE_Synchro.time,SPIKY_loop_results_str{k}.SPIKE_Synchro.sprofile,'color',cm(2,:));
end

figure;hold on;
tmp_mean = mean(SPIKY_interp_rnd);
tmp_se = std(SPIKY_interp_rnd)/sqrt(length(SPIKY_loop_results_rnd));
inan = isnan(tmp_mean) & isnan(tmp_se);tmp_mean(inan) = 0;tmp_se(inan) = 0;

x = [irange(1),irange,irange(end)];
yy1 = [tmp_mean(1),tmp_mean+tmp_se,tmp_mean(end)];
yy2 = [tmp_mean(1),tmp_mean-tmp_se,tmp_mean(end)];
x = [x,x(end:-1:1)];        % repeat x values
yy = [yy1,yy2(end:-1:1)];   % vector of upper & lower boundaries
fill(x,yy,'b','facecolor',cm(1,:),'edgecolor','none')    % fill area defined by x & yy in blue
hr = plot(irange,tmp_mean,'color',cm(1,:));

tmp_mean = mean(SPIKY_interp_str);
tmp_se = std(SPIKY_interp_str)/sqrt(length(SPIKY_loop_results_str));
inan = isnan(tmp_mean) & isnan(tmp_se);tmp_mean(inan) = 0;tmp_se(inan) = 0;

x = [irange(1),irange,irange(end)];
yy1 = [tmp_mean(1),tmp_mean+tmp_se,tmp_mean(end)];
yy2 = [tmp_mean(1),tmp_mean-tmp_se,tmp_mean(end)];
x = [x,x(end:-1:1)];        % repeat x values
yy = [yy1,yy2(end:-1:1)];   % vector of upper & lower boundaries
fill(x,yy,'b','facecolor',cm(2,:),'edgecolor','none')    % fill area defined by x & yy in blue
hs = plot(irange,tmp_mean,'color',cm(2,:));

legend([hr,hs],{'Random','Structured'});xlabel('Time (ms)');ylabel('SPIKE Synchronous metric');title('Cluster synchronicity ');


%% Differentiate UP states:
bgvth = {};
bgvtsh = {};
for ru=1
    for c=1:700
        bgvth{c,ru} = histcounts(batch_bg{c,ru},edges);
        bgvtsh{c,ru} = histcounts(smooth(batch_bg{c,ru},100,'moving'),edges);
    end
end
%%
figure(1);hold on;
figure(2);hold on;
for c=1:100
    figure(1);
    plot(edges(1:end-1),bgvth{c,ru});
    plot(edges(1:end-1),bgvtsh{c,ru},'g');
    f = fit(edges(1:end-1)',bgvtsh{c,ru}','gauss2');
    plot(f,edges(1:end-1)',bgvtsh{c,ru}');
    plot([f.b1-(f.c1*3),f.b1-(f.c1*3)],[0,max(max(bgvtsh{c,ru}))],'b');
    plot([f.b1+(f.c1*3),f.b1+(f.c1*3)],[0,max(max(bgvtsh{c,ru}))],'b');
    plot([f.b2-(f.c2*3),f.b2-(f.c2*3)],[0,max(max(bgvtsh{c,ru}))],'k');
    plot([f.b2+(f.c2*3),f.b2+(f.c2*3)],[0,max(max(bgvtsh{c,ru}))],'k');
    figure(2);
    plot(batch_bg{c,ru});
    plot([0,20000],[f.b1-(f.c1*3),f.b1-(f.c1*3)],'b');
    plot([0,20000],[f.b1+(f.c1*3),f.b1+(f.c1*3)],'b');
    plot([0,20000],[f.b2-(f.c2*3),f.b2-(f.c2*3)],'k');
    plot([0,20000],[f.b2+(f.c2*3),f.b2+(f.c2*3)],'k');
    plot([0,20000],[-66,-66],'r');
    axis([0, 20000, -80, 60]);
    pause();figure(1);cla;figure(2);cla;
end
title('Structured');


% figure(1);hold on;
% figure(2);hold on;
% for c=1:700
%     figure(1);
%     plot(edges(1:end-1),bgvth{c,ru});
%     plot(edges(1:end-1),bgvtsh{c,ru},'r');
%     figure(2);
%     plot(batch_bg{c,ru});
%     axis([0, 20000, -80, 60]);
%     pause();figure(1);cla;figure(2);cla;
% end
% title('Structured');
%% Calculate Voltage trace distributions:
% These are correct after all!!

f_delay_str = cell(700,run.nruns);
gof_delay_str = cell(700,run.nruns);
histo_str = cell(700,run.nruns);
edges_str = cell(700,run.nruns);

f_delay_rnd = cell(700,run.nruns);
gof_delay_rnd = cell(700,run.nruns);
histo_rnd = cell(700,run.nruns);
edges_rnd = cell(700,run.nruns);

edgesint = 0.5;
for ru=1:run.nruns
    [ru]
    parfor c=1:700
%         % For Delay period:
        data = smooth(batch_rnd{c,ru}(1500:end),100,'moving');
        minr = floor(min(data));
        maxr = ceil(max(data));
        n = length(data);
        edges_rnd{c,ru} = minr:edgesint:maxr;
        histo_rnd{c,ru} = histc(data,edges_rnd{c,ru})/(n*edgesint);histo_rnd{c,ru} = histo_rnd{c,ru}';

        [f_delay_rnd{c,ru}, gof_delay_rnd{c,ru}] = fit(edges_rnd{c,ru}',histo_rnd{c,ru}','gauss2');
    end
end
for ru=1:run.nruns
    [ru]
    parfor c=1:700
%         % For Delay period:
        data = smooth(batch_str{c,ru}(1500:end),100,'moving');
        minr = floor(min(data));
        maxr = ceil(max(data));
        n = length(data);
        edges_str{c,ru} = minr:edgesint:maxr;
        histo_str{c,ru} = histc(data,edges_str{c,ru})/(n*edgesint);histo_str{c,ru} = histo_str{c,ru}';

        [f_delay_str{c,ru}, gof_delay_str{c,ru}] = fit(edges_str{c,ru}',histo_str{c,ru}','gauss2');
    end
end


%% Cluster Vm to detect UP states.

MU_delay_rnd = zeros(700,run.nruns);
SIGMA_delay_rnd = zeros(700,run.nruns);
MU_delay_str = zeros(700,run.nruns);
SIGMA_delay_str = zeros(700,run.nruns);

for ru=1:run.nruns
    for c=1:run.nPC
        % Delay period:
        MU_1 = f_delay_rnd{c,ru}.b1;
        MU_2 = f_delay_rnd{c,ru}.b2;
        SIGMA_1 = f_delay_rnd{c,ru}.c1 / sqrt(2);
        SIGMA_2 = f_delay_rnd{c,ru}.c2 / sqrt(2);
        if MU_1 < MU_2
            MU_delay_rnd(c,ru) = MU_1;
            SIGMA_delay_rnd(c,ru) = SIGMA_1;
        else
            MU_delay_rnd(c,ru) = MU_2;
            SIGMA_delay_rnd(c,ru) = SIGMA_2;
        end
        MU_1 = f_delay_str{c,ru}.b1;
        MU_2 = f_delay_str{c,ru}.b2;
        SIGMA_1 = f_delay_str{c,ru}.c1 / sqrt(2);
        SIGMA_2 = f_delay_str{c,ru}.c2 / sqrt(2);
        if MU_1 < MU_2
            MU_delay_str(c,ru) = MU_1;
            SIGMA_delay_str(c,ru) = SIGMA_1;
        else
            MU_delay_str(c,ru) = MU_2;
            SIGMA_delay_str(c,ru) = SIGMA_2;
        end
    end
end

rmse_rnd = zeros(700,run.nruns);
rmse_str = zeros(700,run.nruns);
for ru=1:run.nruns
    for c=1:run.nPC
        rmse_rnd(c,ru) = gof_delay_rnd{c,ru}.rmse;
        rmse_str(c,ru) = gof_delay_str{c,ru}.rmse;
    end
end

% save(sprintf('X:\\Documents\\Glia\\dataParsed2Matlab\\updatedStimGABAb01NEWBGST_mVfits20_SN%d_matrices.mat',run.sn),'MU_delay_rnd','SIGMA_delay_rnd','MU_delay_str','SIGMA_delay_str','rmse_rnd','rmse_str','-v7.3');

% % ONLY FOR THIS RUN:
% MU_delay_rnd = MU_delay_rnd(:,1:66);
% SIGMA_delay_rnd = SIGMA_delay_rnd(:,1:66);
% rmse_rnd = rmse_rnd(:,1:66);

% Plot GoF:
figure;hold on;
mr = 0:0.001:0.02;
plot(mr(2:end),histcounts(rmse_rnd,mr)/numel(rmse_rnd))
plot(mr(2:end),histcounts(rmse_str,mr)/numel(rmse_str),'r')
% I need a better criterion:
rm_GoF_rnd = rmse_rnd(:) > 0.014 ;
rm_GoF_str = rmse_str(:) > 0.014 ;


% Cluster by mixture of gaussians:
nG = 4; % Number of gaussians
X_rnd =[reshape(SIGMA_delay_rnd,1,[])',reshape(MU_delay_rnd,1,[])'];
% % remove spurious points (fitting errors; to check them either way): 
% rm_MU = X_rnd(:,2) < -100 ;
% rm_SI = X_rnd(:,1) > 6 ;
% X_rnd( rm_MU | rm_SI ,: )  = [];
% Can I also do it with GoF criterion?
X_rnd(rm_GoF_rnd,:) = [];

gmfit_rnd = fitgmdist(X_rnd,nG);
clusterX_rnd = cluster(gmfit_rnd,X_rnd);
figure;hold on;
gscatter(X_rnd(:,1),X_rnd(:,2),clusterX_rnd);
scatter(gmfit_rnd.mu(:,1),gmfit_rnd.mu(:,2),'ko');
axis([0, 4, -80, -40]);
title(sprintf('Random:SN%d (All cells/runs)',run.sn));
xlabel('G1 ?');ylabel('G1 ?');

X_str =[reshape(SIGMA_delay_str,1,[])',reshape(MU_delay_str,1,[])'];
% % remove spurious points (fitting errors; to check them either way): 
% rm_MU = X_str(:,2) < -100 ;
% rm_SI = X_str(:,1) > 6 ;
% X_str( rm_MU | rm_SI ,: )  = [];
% Can I also do it with GoF criterion?
X_str(rm_GoF_str,:) = [];
gmfit_str = fitgmdist(X_str,nG);
clusterX_str = cluster(gmfit_str,X_str);
figure;hold on;
gscatter(X_str(:,1),X_str(:,2),clusterX_str);
scatter(gmfit_str.mu(:,1),gmfit_str.mu(:,2),'ko');
axis([0, 4, -80, -40]);
title(sprintf('Structured:SN%d (All cells/runs)',run.sn));
xlabel('G1 ?');ylabel('G1 ?');

%% Visualize (inspect) fitting:
% UPNum_str = zeros(100,1);
% UPFF_str = cell(100,1);
% UPDur_str = cell(100,1);
% UPNum_rnd = zeros(100,1);
% UPFF_rnd = cell(100,1);
% UPDur_rnd = cell(100,1);

%Vars for plotting:
p='plot';
n=18502;
edgesint = 0.1;
edges = -80:edgesint:-30;
if (strcmp(p,'plot'))
    figure(1);hold on;
    figure(2);hold on;
end

for ru=4
    for c=sc_rnd'
        [ru,c]
%         if f_delay_str{c,ru}.b1 < f_delay_str{c,ru}.b2
%             rmp = f_delay_str{c,ru}.b1;
%             ut = f_delay_str{c,ru}.c1*3;
%             lt = 2;
%         else
%             rmp = f_delay_str{c,ru}.b2;
%             ut = f_delay_str{c,ru}.c2*3;
%             lt = 2;
%         end
        if f_delay_rnd{c,ru}.b1 < f_delay_rnd{c,ru}.b2
            rmp = f_delay_rnd{c,ru}.b1;
            ut = f_delay_rnd{c,ru}.c1*3;
            lt = 2;
        else
            rmp = f_delay_rnd{c,ru}.b2;
            ut = f_delay_rnd{c,ru}.c2*3;
            lt = 2;
        end
        

%         if f_delay_str{c,ru}{1}.mu1 < f_delay_str{c,ru}{1}.mu2
%             rmp = f_delay_str{c,ru}{1}.mu1;
%             ut = f_delay_str{c,ru}{1}.sigma1*3;
%             lt = 2;
%         else
%             rmp = f_delay_str{c,ru}{1}.mu2;
%             ut = f_delay_str{c,ru}{1}.sigma2*3;
%             lt = 2;
%         end

%         if f_delay_rnd{c,ru}{1}.mu1 < f_delay_rnd{c,ru}{1}.mu2
%             rmp = f_delay_rnd{c,ru}{1}.mu1;
%             ut = f_delay_rnd{c,ru}{1}.sigma1*3;
%             lt = 2;
%         else
%             rmp = f_delay_rnd{c,ru}{1}.mu2;
%             ut = f_delay_rnd{c,ru}{1}.sigma2*3;
%             lt = 2;
%         end
        
        if (strcmp(p,'plot')) % plot
            find_up_states(tmp{c,ru}(:), lt, ut, rmp, 1, p, gca);
            axis([0, 18500, -80, 60]);
%             tmp = histcounts(smooth(batch_str{c,ru}(:),100,'moving'),edges)/(n*edgesint);
%             
%             fitdata = fcell{1}(f_delay_str{c,ru}{1}.a1,f_delay_str{c,ru}{1}.mu1,f_delay_str{c,ru}{1}.sigma1,f_delay_str{c,ru}{1}.a2,f_delay_str{c,ru}{1}.mu2,f_delay_str{c,ru}{1}.sigma2,edges);
%             
%             figure(1);
%             scatter(edges(1:end-1),tmp,'.k');plot(edges,fitdata,'b');
%             legend({'Original',sprintf('Gauss fit (%f)',gof_delay_str{c,ru}{1}.rmse)});
%             plot([f_delay_str{c,ru}{1}.mu1+f_delay_str{c,ru}{1}.sigma1*3,f_delay_str{c,ru}{1}.mu1+f_delay_str{c,ru}{1}.sigma1*3],[0,max(tmp)],'g');
%             plot([f_delay_str{c,ru}{1}.mu1-f_delay_str{c,ru}{1}.sigma1*3,f_delay_str{c,ru}{1}.mu1-f_delay_str{c,ru}{1}.sigma1*3],[0,max(tmp)],'g');
%             plot([f_delay_str{c,ru}{1}.mu2+f_delay_str{c,ru}{1}.sigma2*3,f_delay_str{c,ru}{1}.mu2+f_delay_str{c,ru}{1}.sigma2*3],[0,max(tmp)],'g');
%             plot([f_delay_str{c,ru}{1}.mu2-f_delay_str{c,ru}{1}.sigma2*3,f_delay_str{c,ru}{1}.mu2-f_delay_str{c,ru}{1}.sigma2*3],[0,max(tmp)],'g');
%             plot(edges,normpdf(edges,f_delay_str{c,ru}{1}.mu1,f_delay_str{c,ru}{1}.sigma1),'c');
%             plot(edges,normpdf(edges,f_delay_str{c,ru}{1}.mu2,f_delay_str{c,ru}{1}.sigma2),'c');
            set(gcf,'Position',[  257   558   983   420]);
            pause();figure(1);cla;figure(2);cla;
        else
%             [~,UPDur_str{c,ru},UPs_str] = find_up_states(batch_str{c,ru}(1500:end), lt, ut, rmp, 1);
%             [~,UPDur_rnd{c,ru},UPs_rnd] = find_up_states(batch_rnd{c,ru}(1500:end), lt, ut, rmp, 1);
        end % plot

%         UPNum_str(c,ru) = length(UPDur_str{c,ru});
%         UPNum_rnd(c,ru) = length(UPDur_rnd{c,ru});
%         if ~isempty(UPs_str)
%             UPFF_rnd{c,ru} = (cellfun(@(x) advanced_spike_count(x, -10, 0),UPs_rnd)./(UPDur_rnd{c,ru}*0.001)')';
%             UPFF_str{c,ru} = (cellfun(@(x) advanced_spike_count(x, -10, 0),UPs_str)./(UPDur_str{c,ru}*0.001)')';
%         end
    end
end

%% FF of spikes uppon Up states
edges = 0:5:80;
%Remove zero frequency; explodes graph:
data = cell2mat(reshape(UPFF_rnd,1,[]));
histo = histcounts(data,edges)/length(data);
figure;hr=plot(edges(2:end-1),histo(2:end));
data = cell2mat(reshape(UPFF_str,1,[]));
histo = histcounts(data,edges)/length(data);
hold on;hs=plot(edges(2:end-1),histo(2:end),'r');
legend([hr,hs],{'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Relative Frequency');title('Firing rate uppon UP state');
%Duration of Up states
edges = 20:20:1000;
data = cell2mat(reshape(UPDur_rnd,1,[]));
histo = histcounts(data,edges)/length(data);
figure;hr=plot(edges(1:end-1),histo);
data = cell2mat(reshape(UPDur_str,1,[]));
histo = histcounts(data,edges)/length(data);
hold on;hs=plot(edges(1:end-1),histo,'r');
legend([hr,hs],{'Random','Structured'});xlabel('Duration in miliseconds');ylabel('Relative Frequency');title('Duration of UP state');
% Number of UP states:
edges = 0:1:50;
tmp_rnd = UPNum_rnd(UPNum_rnd~=0);
tmp_str = UPNum_str(UPNum_str~=0);
tmp = [histcounts(tmp_rnd(:),edges)'/length(tmp_rnd),histcounts(tmp_str(:),edges)'/length(tmp_str)];
cm = get(groot,'DefaultAxesColorOrder');
% figure;hb = bar(edges(2:end-1),tmp(2:end,:));
figure;hb = plot(edges(2:end-1),tmp(2:end,:));
% hb(1).FaceColor = cm(1,:);
% hb(2).FaceColor = cm(2,:);
legend({'Random','Structured'},'Location','NorthWest');xlabel('Number of UP states');ylabel('Relative Frequency');title('Number of UP states');


%%

figure(1);hold on;
figure(2);hold on;
for ru = 1:1
    for c=1:length(sc_str)
        [c,ru]
        figure(1);
        plot(edges(1:end-1),vth{c,ru});
        plot(edges(1:end-1),vtsh{c,ru},'g');
        if (0)% gauss2
        f = fit(edges(1:end-1)',vtsh{c,ru}','gauss2');
        if f.b1 < f.b2
            rmp = f.b1;
            ut = (f.b1+(f.c1*3))-rmp;
            lt = 2;
        else
            rmp = f.b2;
            ut = (f.b2+(f.c2*3))-rmp;
            lt = 2;
        end
        plot(f,edges(1:end-1)',vtsh{c,ru}');
        plot([f.b1-(f.c1*3),f.b1-(f.c1*3)],[0,max(max(vtsh{c,ru}))],'b');
        plot([f.b1+(f.c1*3),f.b1+(f.c1*3)],[0,max(max(vtsh{c,ru}))],'b');
        plot([f.b1,f.b1],[0,max(max(vtsh{c,ru}))],'b');
        plot([f.b2-(f.c2*3),f.b2-(f.c2*3)],[0,max(max(vtsh{c,ru}))],'k');
        plot([f.b2+(f.c2*3),f.b2+(f.c2*3)],[0,max(max(vtsh{c,ru}))],'k');
        plot([f.b2,f.b2],[0,max(max(vtsh{c,ru}))],'k');
        elseif (1) % gauss
            f = fit(edges(1:end-1)',vtsh{c,ru}','a1*exp(-((x-b1)/c1)^2)');
            rmp = f.b1;
            ut = (f.b1+(f.c1*3))-rmp;
            lt = 2;
            plot(f,edges(1:end-1)',vtsh{c,ru}');
            plot([f.b1-(f.c1*3),f.b1-(f.c1*3)],[0,max(max(vtsh{c,ru}))],'b');
            plot([f.b1+(f.c1*3),f.b1+(f.c1*3)],[0,max(max(vtsh{c,ru}))],'b');
        elseif (0) %skew normal
        end
        figure(2);
%         UT(c,ru) = ut;
%         LT(c,ru) = lt;
%         RMP(c,ru) = rmp;
        [S,D,UPs] = find_up_states(batch_str{sc_str(c),ru}(1500:end), lt, ut, rmp, 0.2, 'plot',gca );
        axis([0, 18500, -80, 60]);
        pause();figure(1);cla;figure(2);cla;
    end
end

% figure;
% scatter(UT(:),RMP(:));
% xlabel('Upper Treshold');ylabel('Resting Membrane Potential');

X = [UT(:),RMP(:)];
gmfit = fitgmdist(X,3);
clusterX = cluster(gmfit,X);
figure;hold on;
h = gscatter(X(:,1),X(:,2),clusterX);

% ezcontour(@(x1,x2)pdf(gmfit,[x1 x2]),get(gca,{'XLim','YLim'}))
% title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
% legend(h,'Model 0','Model1')
% hold off

%% Distribution example:
% What's the difference?


minr = 10;
maxr = 80;
rng default;  % For reproducibility
n = 10000;
histobinint = 0.2;
databinint = 0.5;

dataX=minr:databinint:maxr;
histoX=minr:histobinint:maxr;
dataseed = rand(1,n);
pdf = mygauss2(0.5,50,2,0.5,45,3,dataX);  % Simulate heights
cdf = cumsum(pdf);cdf = cdf/cdf(end);[cdf,mask]=unique(cdf);
data = interp1(cdf(mask), dataX(mask), dataseed);
figure;hold on;
% scatter(data,zeros(1,n),'.');
plot(dataX,pdf,'b')
% [mu,s,muci,sci] = normfit(height)

histo = histcounts(data,histoX)/(n*histobinint);

scatter(histoX(1:end-1),histo,'g');

myfittype = fittype(mygauss2)
% myfittype = fittype('(1/(s*sqrt(2*pi)))*exp(-((x-mu)/(2*s)).^2)',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'mu','s'});
myfitoptions = fitoptions(myfittype);
myfitoptions.Robust = 'on';
myfitoptiond.Algorithm = 'Levenberg-Marquardt';
myfitoptions.Lower = [0.01, 20, 0.01, 0.01, 20, 0.01];
myfitoptions.Upper = [1, 60, 6, 1, 60, 6];
myfitoptions.StartPoint = [0.5, 30, 1, 0.5, 50, 1];
myfitoptions.Display = 'iter';
myfitoptions.TolFun = 1.0e-09;
[f, gof] = fit(histoX(1:end-1)',histo',myfittype,myfitoptions);

plot(f,'k');

% fprintf('In sigma 1 %.3f (%.3f)\n',68.2,sum((height<52)&(height>48))/n)
% fprintf('In sigma 2 %.3f (%.3f)\n',95.4,sum((height<54)&(height>46))/n)
% fprintf('In sigma 2 %.3f (%.3f)\n',99.6,sum((height<56)&(height>44))/n)


%%
f = fit([1:length(anss)]',anss','gauss2')
figure;plot(f,[1:length(anss)]',anss')

% figure;hold on;
% for c=1:100
% plot(edges(1:end-1),vth{c,ru});
% plot(edges(1:end-1),vtsh{c,ru},'r');
% % scatter(vth{c,2},1,'+');
% pause();cla;
% end
% title('Structured');
%%
plot_ifr(batch(run.StimMat_str(:)),tstop);

%%
figure;hold on;
for k=1:700
    spikes = batch{run.StimMat_str(k)};
    if ~isempty(spikes)
        scatter(spikes,ones(1,length(spikes))*k,'+');
    end
end