close all;clear all;clc;
% Load run properties:
% load('multi_core_R/CRUN.mat');
load('../../../CRUN.mat');
clear PCcells mycell;
Runs = 1;

for r = Runs
    for c=1:nPC
        mycell = ncell(load(sprintf('multi_core/%d_%d.txt',c-1,r-1)),10);
        mycell.clusterID = labels(c,Sid);
%         mycell=mycell.binSpikes(50);
%        
        PCcells{c,r}=mycell.hasPersistent(1000,8,4000);
        
        %     plot( PCcells{i}.mv,'linewidth',1,'color',rand(1,3));
        %     pause;
        %     cla
    end
end

imagesc([cellfun(@(x) x.persistentActivity,PCcells),StimVect])

% for r = Runs
%     for c=1:nPC
%         plot( PCcells{c,r}.mv,'linewidth',1,'color',rand(1,3));hold on;
%         axis([0,length(PCcells{c,r}.mv),-70,70 ]);hold on;
%         pause;
%         cla;
%     end
% end


for r = 4
    for cl=1:NC(Sid) 
        figure(); hold on;
        set(gcf,'numbertitle','off','name', sprintf('Cluster #%d',cl) );
        for c=1:nPC
            PCcells{c,r}.clusterID
            if(PCcells{c,r}.clusterID == cl)
                plot( PCcells{c,r}.mv,'linewidth',1,'color',rand(1,3));
            end
        end
        axis([0,length(PCcells{c,r}.mv),-70,70 ]);hold on;
    end
end

%%
figure(); hold on;
axis([0,length(parallel),-70,70 ])
%%
close all;clear all;clc;
for i=[1,24:48]
    parallel = load(sprintf('multi_core/%d_0.txt',i-1));
    single = load(sprintf('single_core/PCsomaV_%d_run_1.txt',i-1));
    plot(parallel,'linewidth',1);figure(gcf);
    hold on;plot(single,'r');
%     axis([0,length(parallel),-70,70 ])
%     all((parallel-single)<1e-4)
%     mean(parallel-single)
    pause;
    cla;
end
%%
for i=0:6
    parallel = load(sprintf('multi_core/%d_0.txt',40+i));
    single = load(sprintf('single_core/PVsomaV_%d_run_1.txt',i));
    plot(parallel,'linewidth',6);figure(gcf);
    hold on;plot(single,'r');
%     mean(parallel-single)
    pause;
    cla;
end
%%
close all;clear all;clc;
for i=1:40
    parallel = load(sprintf('multi_core/%d_0.txt',i-1));
    plot(parallel,'linewidth',6);figure(gcf);
    pause;
    cla;
end