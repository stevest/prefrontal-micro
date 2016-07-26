function getMeanStateFrequencyPerQ(run,p)
%% Compare states across configurations:
% VARPID = 25;
x_str = []; g_str = [];
x_rnd = []; g_rnd = [];
for Qi = 1:length(p.Qseq)
    p.Qi = Qi;
    
    % Stats for Random configuration:
    p.activeconfig = 'rnd';
    eval( sprintf('stats_%s = getConfigurationStats(run,p);',p.activeconfig) );
    
    % Stats for Structured configuration:
    p.activeconfig = 'str';
    eval( sprintf('stats_%s = getConfigurationStats(run,p);',p.activeconfig) );    
    
    % EXCLUDE ZERO STATE: (is outlier in boxplots)
    x_rnd = [x_rnd ; reshape(stats_rnd.maxfreqstates(2:end),[],1) ] ;
    g_rnd = [g_rnd ; ones(length(stats_rnd.maxfreqstates)-1,1)*Qi ] ;
    x_str = [x_str ; reshape(stats_str.maxfreqstates(2:end),[],1) ] ;
    g_str = [g_str ; ones(length(stats_str.maxfreqstates)-1,1)*Qi ] ;
    
end
figure;hold on;
boxplot(x_rnd,g_rnd);
scatter(g_rnd,x_rnd,10,'k','filled');
set(gca,'XtickLabel',p.Qseq);
title(sprintf('Random Q=All (SN%d)',run.sn));
xlabel('Different Qs');ylabel('Frequency of States (pooled)');


figure;hold on;
boxplot(x_str,g_str);
scatter(g_str,x_str,10,'k','filled');
set(gca,'XtickLabel',p.Qseq);
title(sprintf('Structured Q=All (SN%d)',run.sn));
xlabel('Different Qs');ylabel('Frequency of States (pooled)');
end
