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
sc_rnd = run.StimMat_rnd(:,1);
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

%% pyramidal CV ONLY FOR STIMULATED CLUSTER AND DELAY PERIOD!
FF_stim_rnd = zeros(100,80);
FF_delay_rnd = zeros(100,80);
CV_stim_rnd = zeros(100,80);
CV_delay_rnd = zeros(100,80);
FF_stim_str = zeros(100,80);
FF_delay_str = zeros(100,80);
CV_stim_str = zeros(100,80);
CV_delay_str = zeros(100,80);
for ru=1:72
    for c=1:100
%             [~,st_rnd] = advanced_spike_count(batch_rnd{sc_rnd(c),ru}(1500:end), -10, 0);
        st_rnd = batch_rnd_spikes{sc_rnd(c),ru};
        st_stim_rnd = st_rnd(st_rnd<=1500);
        st_delay_rnd = st_rnd(st_rnd>1500);
        ISI_stim_rnd = diff(st_stim_rnd);
        ISI_delay_rnd = diff(st_delay_rnd);
        FF_stim_rnd(c,ru) = length(st_stim_rnd)/18.5;
        FF_delay_rnd(c,ru) = length(st_delay_rnd)/18.5;
        CV_stim_rnd(c,ru) = std(ISI_stim_rnd)/mean(ISI_stim_rnd);
        CV_delay_rnd(c,ru) = std(ISI_delay_rnd)/mean(ISI_delay_rnd);
%         [~,st_str] = advanced_spike_count(batch_str{sc_str(c),ru}(1500:end), -10, 0);
        st_str = batch_str_spikes{sc_str(c),ru};
        st_stim_str = st_str(st_str<=1500);
        st_delay_str = st_str(st_str>1500);
        ISI_stim_str= diff(st_stim_str);
        ISI_delay_str = diff(st_delay_str);
        FF_stim_str(c,ru) = length(st_stim_str)/18.5;
        FF_delay_str(c,ru) = length(st_delay_str)/18.5;
        CV_stim_str(c,ru) = std(ISI_stim_str)/mean(ISI_stim_str);
        CV_delay_str(c,ru) = std(ISI_delay_str)/mean(ISI_delay_str);
    end
end

% CV for stimulated cluster (delay period):
edges = 0:0.1:2;
tmp_delay = [histcounts(CV_delay_rnd(:),edges)',histcounts(CV_delay_str(:),edges)']/80;
figure;hold on;plot(edges(2:end-1),tmp_delay(2:end,:));
cm = get(groot,'DefaultAxesColorOrder');
tmp_stim = [histcounts(CV_stim_rnd(:),edges)',histcounts(CV_stim_str(:),edges)']/80;
plot(edges(2:end-1),tmp_stim(2:end,1),'color',cm(1,:));
plot(edges(2:end-1),tmp_stim(2:end,2),'color',cm(2,:));

dmax=max([tmp_delay(:);tmp_stim(:)]);dmax=dmax+dmax*0.1;
scatter(nanmean([CV_delay_rnd(:),CV_delay_str(:)]),[1,1]*dmax,[],cm(1:2,:),'v','fill');
scatter(nanmean([CV_stim_rnd(:),CV_stim_str(:)]),[1,1]*dmax,[],cm(1:2,:),'v','fill');

legend({'Random','Structured'});xlabel('Coefficient of Variation (CV)');ylabel('Count (#)');title('Histogram of CV (Cluster in delay period)');

% FF for stimulated cluster (delay period):
edges = 0:2:60;
tmp = [histcounts(FF_delay_rnd(:),edges)',histcounts(FF_delay_str(:),edges)']/80;
figure;plot(edges(2:end-1),tmp(2:end,:));
legend({'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Count (#)');title('Histogram of FF (Cluster in delay period)');
set(gca,'XScale','log');

%% Compute Synchronicity:
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

%% Recruitement:
FFr_rnd = zeros(100,72);
FFr_str = zeros(100,72);

for ru=1:72
    for c=1:600
        st_rnd = batch_rnd_spikes{nsc_rnd(c),ru};
        st_rnd = st_rnd(st_rnd>1500);
        FFr_rnd(c,ru) = length(st_rnd)/18.5;

        st_str = batch_str_spikes{nsc_str(c),ru};
        st_str = st_str(st_str>1500);
        FFr_str(c,ru) = length(st_str)/18.5;
    end
end


edges = 0:0.2:8;
tmp = [histcounts(FFr_rnd(FFr_rnd~=0),edges)',histcounts(FFr_str(FFr_str~=0),edges)']/(72*600);
% figure;plot(edges,histcounts(FFr_rnd(FFr_rnd~=0),edges))
% figure;plot(edges,histcounts(FFr_str(FFr_str~=0),edges))
figure;hold on;plot(edges(1:end-1),tmp)
legend({'Random','Structured'});xlabel('Firing Frequency (Hz)');ylabel('Relative Frequency');title('Non-Stimulated Cells Firing Rate');
set(gca,'XScale','log');
plot([0.4,0.4],[0,0.04])

% Sort by mean activity (per cell):
[MEAN_rnd,IDX_rnd] = sort(mean(FFr_rnd,2),'descend');
[MEAN_str,IDX_str] = sort(mean(FFr_str(:,[1:5,7:19,21:73,75:80]),2),'descend');

% Set limit in background activity (1Hz):
figure;hold on;
plot(MEAN_rnd(1:40));
plot(MEAN_str(1:40),'r');
xlabel('Cell ID (#)');ylabel('Mean Firing Frequency');

cmmax = max(max([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,[1:5,7:19,21:73,75:80])]));
cmmin = min(min([FFr_rnd(IDX_rnd,:),FFr_str(IDX_str,[1:5,7:19,21:73,75:80])]));

figure;imagesc(FFr_rnd(IDX_rnd,:));
caxis manual
caxis([cmmin cmmax]);
cm = hot(1000);
cm(1,:) = [0,0,0];
colormap(cm);title('Random');
figure;imagesc(FFr_str(IDX_str,[1:5,7:19,21:73,75:80]));
caxis manual
caxis([cmmin cmmax]);
colormap(cm);title('Structured');

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
% sc_str = stimdata_str(:,1);
% sc_rnd = stimdata_rnd(:,1);
sc_str = run.StimMat_str(:,1);
nsc_str = find(~ismember(1:700,run.StimMat_str(:,1)));
sc_rnd = run.StimMat_rnd(:,1);
nsc_rnd = find(~ismember(1:700,run.StimMat_rnd(:,1)));

% n=18502;
% edgesint = 0.5;
% edges = -80:edgesint:-40;
% vth = {};
% vtsh = {};

% fcell{1} = @(mu,sigma,x) exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma); %gauss
% fcell{2} = @(a1,mu1,sigma1,a2,mu2,sigma2,x) a1*fcell{1}(mu1,sigma1,x) + a2*fcell{1}(mu2,sigma2,x); %gauss2
% fcell{3} = @(mu,a,sigma,x) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,mu,sigma);
clear fcell fparams;
mygauss = @(mu,sigma,x) exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma); %gauss
fcell{1} = @(a1,mu1,sigma1,a2,mu2,sigma2,x) a1*mygauss(mu1,sigma1,x) + a2*mygauss(mu2,sigma2,x); %gauss2

fparams{1,1} = [0.01, -80, 0.001, 0.01, -80, 0.001];
fparams{1,2} = [0.99, -30, 10, 0.99, -30, 10];
% Change fparams to give min max bounds to randomly select startpoint:
% fparams{1,1} = [-66, 2];
% fparams{1,2} = [-80, 0.001];
% fparams{1,3} = [-30, 10];
% fparams{2,1} = [0.5, -76, 2, 0.5, -36, 2];
% fparams{2,2} = [0.01, -80, 0.001, 0.01, -80, 0.001];
% fparams{2,3} = [0.99, -30, 10, 0.99, -30, 10];
% fparams{3,1} = [-66, 0, 2];
% fparams{3,2} = [-80, -20, 0.001];
% fparams{3,3} = [-30, 20, 10];

% fparams{1,1} = [-80, 0.001];
% fparams{1,2} = [-30, 10];
% fparams{2,1} = [0.01, -80, 0.001, 0.01, -80, 0.001];
% fparams{2,2} = [0.99, -30, 10, 0.99, -30, 10];
% fparams{3,1} = [-80, -20, 0.001];
% fparams{3,2} = [-30, 20, 10];

parpool(8);



f_delay_str = cell(700,72);
gof_delay_str = cell(700,72);
% histo_str = cell(700,72);
% edges_str = cell(700,72);

% f_delay_rnd = cell(700,72);
% gof_delay_rnd = cell(700,72);
% histo_rnd = cell(700,72);
% edges_rnd = cell(700,72);

edgesint = 0.5;
for ru=1:72
    [ru]
    parfor c=1:700%sc_rnd'
%         [ru,c]
%         vth{c,ru} = histcounts(batch_str{sc_str(c),ru}(1500:end),edges)/(n*edgesint);
%         vtsh{c,ru} = histcounts(smooth(batch_str{sc_str(c),ru}(1500:end),100,'moving'),edges)/(n*edgesint);
%         % For Delay period:
%         data = smooth(batch_str{c,ru}(1500:end),100,'moving');
%         minr = floor(min(data));
%         maxr = ceil(max(data));
%         n = length(data);
%         edges_str{c,ru} = minr:edgesint:maxr;
%         histo_str{c,ru} = histcounts(data,edges_str{c,ru})/(n*edgesint);histo_str{c,ru} = histo_str{c,ru}';
        
%             [f_delay_str{c,ru}, gof_delay_str{c,ru}] = best_fitting_function(smooth(batch_str{c,ru}(1500:end),100,'moving'), fcell, fparams);
%         [f_delay_str{c,ru}, gof_delay_str{c,ru}] = fast_fitting_function(histo_str{c,ru}, fcell, fparams, edges_str{c,ru});
%         [f_delay_rnd{c,ru}, gof_delay_rnd{c,ru}] = fast_fitting_function(histo_rnd{c,ru}, fcell, fparams, edges_rnd{c,ru},10);

%         [f_delay_rnd{c,ru}, gof_delay_rnd{c,ru}] = fit(edges_rnd{c,ru}(1:end-1)',histo_rnd{c,ru},'gauss2');
        [f_delay_str{c,ru}, gof_delay_str{c,ru}] = fit(edges_str{c,ru}(1:end-1)',histo_str{c,ru},'gauss2');

%         % For Stimulation period:
%         [f_stim{c,ru}, gof_stim{c,ru}] = best_fitting_function(smooth(batch_str{c,ru}(500:1500),100,'moving'), fcell, fparams);

%         % For Baseline period:
%         [f_base{c,ru}, gof_base{c,ru}] = best_fitting_function(smooth(batch_bg{c,ru},100,'moving'), fcell, fparams);
    end
end


MU_delay_rnd = zeros(700,72);
SIGMA_delay_rnd = zeros(700,72);
MU_delay_str = zeros(700,72);
SIGMA_delay_str = zeros(700,72);
% MU_base = zeros(700,1);
% SIGMA_base = zeros(700,1);
% MU_stim = zeros(700,1);
% SIGMA_stim = zeros(700,1);
for ru=1:72
    for c=1:700
%         % Baseline period:
%         if f_base{c,ru}{1}.mu1 < f_base{c,ru}{1}.mu2
%             MU_base(c,ru) = f_base{c,ru}{1}.mu1;
%             SIGMA_base(c,ru) = f_base{c,ru}{1}.sigma1;
%         else
%             MU_base(c,ru) = f_base{c,ru}{1}.mu2;
%             SIGMA_base(c,ru) = f_base{c,ru}{1}.sigma2;
%         end

%         % Stimulation period:
%         if f_stim{c,ru}{1}.mu1 < f_stim{c,ru}{1}.mu2
%             MU_stim(c,ru) = f_stim{c,ru}{1}.mu1;
%             SIGMA_stim(c,ru) = f_stim{c,ru}{1}.sigma1;
%         else
%             MU_stim(c,ru) = f_stim{c,ru}{1}.mu2;
%             SIGMA_stim(c,ru) = f_stim{c,ru}{1}.sigma2;
%         end

        % Delay period:
%         if f_delay_rnd{c,ru}{1}.mu1 < f_delay_rnd{c,ru}{1}.mu2
%             MU_delay(c,ru) = f_delay_rnd{c,ru}{1}.mu1;
%             SIGMA_delay(c,ru) = f_delay_rnd{c,ru}{1}.sigma1;
%         else
%             MU_delay(c,ru) = f_delay_rnd{c,ru}{1}.mu2;
%             SIGMA_delay(c,ru) = f_delay_rnd{c,ru}{1}.sigma2;
%         end
        if f_delay_rnd{c,ru}.b1 < f_delay_rnd{c,ru}.b2
            MU_delay_rnd(c,ru) = f_delay_rnd{c,ru}.b1;
            SIGMA_delay_rnd(c,ru) = f_delay_rnd{c,ru}.c1;
        else
            MU_delay_rnd(c,ru) = f_delay_rnd{c,ru}.b2;
            SIGMA_delay_rnd(c,ru) = f_delay_rnd{c,ru}.c2;
        end
        if f_delay_str{c,ru}.b1 < f_delay_str{c,ru}.b2
            MU_delay_str(c,ru) = f_delay_str{c,ru}.b1;
            SIGMA_delay_str(c,ru) = f_delay_str{c,ru}.c1;
        else
            MU_delay_str(c,ru) = f_delay_str{c,ru}.b2;
            SIGMA_delay_str(c,ru) = f_delay_str{c,ru}.c2;
        end
    end
end

sc_str = run.StimMat_str(:,1);
nsc_str = find(~ismember(1:700,run.StimMat_str(:,1)));
sc_rnd = run.StimMat_rnd(:,1);
nsc_rnd = find(~ismember(1:700,run.StimMat_rnd(:,1)));


tmp_MU_IDX_rnd = MU_delay_rnd(:)>-80 ;
tmp_MU_IDX_str = MU_delay_str(:)>-80 ;

figure;
scatter(reshape(SIGMA_delay_rnd(sc_rnd,:),1,[]),reshape(MU_delay_rnd(sc_rnd,:),1,[]),'g.');
axis([0, 4, -80, -40]);
X =[reshape(SIGMA_delay_rnd(sc_rnd,:),1,[])',reshape(MU_delay_rnd(sc_rnd,:),1,[])'];
gmfit_rnd = fitgmdist(X,3);
clusterX_rnd = cluster(gmfit_rnd,X);
figure;hold on;
h = gscatter(X(:,1),X(:,2),clusterX_rnd);
scatter(gmfit_rnd.mu(:,1),gmfit_rnd.mu(:,2),'ko');
axis([0, 4, -80, -40]);

figure;
scatter(reshape(SIGMA_delay_str(sc_str,:),1,[]),reshape(MU_delay_str(sc_str,:),1,[]),'r.');
axis([0, 4, -80, -40]);
% X =[reshape(SIGMA_delay_str(sc_str,:),1,[])',reshape(MU_delay_str(sc_str,:),1,[])'];
% gmfit_str = fitgmdist(X,4);
% clusterX_str = cluster(gmfit_str,X);
figure;hold on;
h = gscatter(X(:,1),X(:,2),clusterX_str);
scatter(gmfit_str.mu(:,1),gmfit_str.mu(:,2),'ko');
axis([0, 4, -80, -40]);


h4=scatter(reshape(SIGMA_delay(sc_str,:),1,[]),reshape(MU_delay(sc_str,:),1,[]),'ko');
h1=scatter(SIGMA_base(:),MU_base(:),'c.');
% h2=scatter(SIGMA_stim(:),MU_stim(:),'g');
legend([h1,h4,h3],{'baseline','stim','delay'});

X =[SIGMA_delay(:),MU_delay(:)];
gmfit = fitgmdist(X,3);
clusterX = cluster(gmfit,X);
figure;hold on;
h = gscatter(X(:,1),X(:,2),clusterX);
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