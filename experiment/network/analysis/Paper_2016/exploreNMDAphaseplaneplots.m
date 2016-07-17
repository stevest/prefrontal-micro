%% load iNMDA and mV for one run:
clear all;close all;
runParams.stc = 4;
runParams.config = 'str';
if strcmp(runParams.config, 'str')
cfgprefix = 'S';
else
    cfgprefix = 'R';
end
runParams.SN = 2;
runParams.runrange = 1:50;
nrun = length(runParams.runrange);
runParams.PID = 25;
runParams.experimentDirStr =     sprintf('SAVENMDAupdatedStimGABAb01NEWBGST_%ss10c%%d_SN%%d_mV_r%%d_PID%%d',cfgprefix);
runParams.experimentFileStr =    [runParams.experimentDirStr,'.mat'];
runParams.args = {runParams.stc-1,runParams.SN,0,runParams.PID};

loadParams.specifics = 'inmda';
loadParams.specificsFileStr =    sprintf('SAVENMDAupdatedStimGABAb01NEWBGST_%ss10c%%d_SN%%d_inmda_r%%d_PID%%d.mat',cfgprefix);

%%
filename = fullfile('\\139.91.162.50\stefanos_synology\','Glia',sprintf(runParams.experimentFileStr, runParams.args{:} ));
load(filename);
nmdafilename = fullfile('\\139.91.162.50\stefanos_synology\','Glia',sprintf(loadParams.specificsFileStr, runParams.args{:} ));
load(nmdafilename);
%%
close all;
figure(1);figure;(2)
for c=29%700
    mV = batch{c,1};
%     mV = nmdatraces(c,:);
    figure(1);plot(mV);
    % mV = interp1(1:length(batch{1,1}),batch{1,1},1:0.1:length(batch{1,1}));
    % figure;plot(1:0.1:length(batch{1,1}),mV,'b');hold on;plot(1:length(batch{1,1}),batch{1,1},'r');
    
%     y = diff(mV);
    y = inmdaBatch{4,c};
    x = mV(1:end-1);
    figure(2);plot(x,y);
%     pause();cla;
end

%%
pathto = '\\139.91.162.90\cluster\stefanos\Documents\GitHub\prefrontal-micro-tmp\experiment\network\DATA\';
% d = dir(pathto);
% files = {d(~[d.isdir]).name};
% 
% nmdatraces = zeros(length(files),25001);
% for fn = 1:length(files)
%    nmdatraces(fn,:) = load([pathto,files{fn}]);
% end
soma = load([pathto,'nmda_loc_70_syn_20_nmdaw_50_ISI_10_soma.txt']);
basal = load([pathto,'nmda_loc_70_syn_20_nmdaw_50_ISI_10_basal.txt']);
binmda = ones(20,500*50+1);
for k=1:20
    binmda(k,:) = load([pathto,sprintf('nmda_loc_70_syn_20_nmdaw_50_ISI_10_basalinmda%d.txt',k-1)]);
end

%% Interactive subplot:
close all;
fig = figure;
c=1;
t = 1;
start = 50*50;
stop = 200*50;
mV = basal;
mV = mV(start:stop);
tstop = length(mV);
subplot(2,2,1);plot(mV);hold on;
vh = scatter(t,mV(t),'filled');
% y = diff(mV);
y = -1*binmda(1,start:stop);
% x = mV(1:end-1);
x = mV;
subplot(2,2,2);plot(x,y);hold on;
ph = scatter(x(t),y(t),'filled');

bl1 = uicontrol('Parent',fig,'Style','text','Position',[50,54,23,23],...
                'String',num2str(t));
b = uicontrol('Parent',fig,'Style','slider','Position',[10,10,1600,10],...
              'value',0, 'min',0, 'max',(tstop-1)/50);
b.Callback = @(obj,ed) clb_phaseplaneplot_update(obj,ed,vh,ph,mV,x,y,bl1);


