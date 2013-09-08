cd('TYPES_Ra_100')
%%
current = 0.05;
figure('Name', '');
for id = 0:55
    for runs=0:8
        load( sprintf('SomaV_%d_%d.txt', id, runs) );
        set(gcf,'Name', sprintf('Morphology #%d',id));
        subplot(3,3,runs+1);hold on;
        title(sprintf( 'IClamp @ %1.2f',current ));
        eval( ['plot( ', sprintf('SomaV_%d_%d', id, runs) , ' );hold on;' ] );
        current = current + 0.05
        
    end
    pause;
    for runs=0:8
        subplot(3,3,runs+1);hold off;cla;
    end
    current =0.05
end
