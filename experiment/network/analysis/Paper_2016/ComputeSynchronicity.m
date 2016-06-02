pathto = '\\139.91.162.90\cluster\stefanos\Documents\Glia\GABAb02_Rs20c0r1';
[batch,tstop] = load_raw_batch(pathto,'spikes');
%%
plot_ifr(batch(run.StimMat_rnd(:)),tstop);

%% Calculate synchronicity:

% construct input array:
spikes = batch(run.StimMat_rnd(:,1));
%Construct parameters struct:
para.tmin = 0;
para.tmax = 20000;
para.dts = 1;
para.num_trains = length(spikes);

m_para.all_measures_string={'SPIKE_Synchro';'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_forward';'PSTH';};  % order of select_measures

para.select_measures      =[1 1 1 0 0 1];  % Select measures (0-calculate,1-do not calculate)

%Compute synchronicity:
SPIKY_loop_results = SPIKY_loop_f_distances(spikes,para) %#ok<NOPTS>

%% Plot!:
plotting=7;           % +1:spikes,+2:dissimilarity profile,+4:dissimilarity matrix
plot_profiles=1;      % 1:all,2:Groups only,3:all+groups

num_trains = para.num_trains;

num_plots=(mod(plotting,2)>0)+(mod(plotting,4)>1)+(mod(plotting,8)>3);
if num_plots>0
    measures=find(para.select_measures);
    for mc=1:length(measures)
        measure=measures(mc);
        measure_var=m_para.all_measures_string{measure};
        measure_name=regexprep(measure_var,'_','-');

        figure(mc); clf
        set(gcf,'Name',measure_name)
        set(gcf,'Units','normalized','Position',[0.0525 0.0342 0.8854 0.8867])
        subplotc=0;

        if mod(plotting,2)>0
            subplotc=subplotc+1;
            subplot(num_plots,1,subplotc)                                      % Spikes
            for trc=1:length(spikes)
                for spc=1:length(spikes{trc})
                    line(spikes{trc}(spc)*ones(1,2),length(spikes)-[trc-1 trc])
                end
            end
            xlim([para.tmin para.tmax])
            ylim([0 length(spikes)])
            if num_trains<10
                set(gca,'YTick',-0.5+(1:num_trains),'YTickLabel',fliplr(1:num_trains))
            else
                set(gca,'YTick',[])
            end
            title ('Spike trains','FontWeight','bold','FontSize',14)
        end

        if mod(plotting,4)>1
            subplotc=subplotc+1;
            subplot(num_plots,1,subplotc)                                      % Dissimilarity profile
            if measure==1                                             % Group profiles for SPIKE-Sync need extra treatment since support is different for each group
                x=SPIKY_loop_results.(measure_var).time;
                y=SPIKY_loop_results.(measure_var).profile;
                hold on
                num_profs=mod(size(y,1),2)+(size(y,1)-mod(size(y,1),2))/2;
                if plot_profiles>1
                    cols='krbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmcrbgmc';
                    for pc=2:num_profs
                        pfindy=find(y(pc+(size(y,1)-mod(size(y,1),2))/2,:));
                        if ~isempty(pfindy)
                            plot(x(pfindy),y(pc,pfindy),cols(pc),'LineWidth',1);
                        end
                    end
                end
                if mod(plot_profiles,2)>0
                    plot(x,y(1,:),'k','LineWidth',1.5)
                end
            elseif measure==2                                              % piecewise constant profiles (ISI) have first to be transformed
                isi_x=SPIKY_loop_results.(measure_var).time;
                isi_y=SPIKY_loop_results.(measure_var).profile;
                plot_y_values=zeros(size(isi_y,1),length(isi_x)*2);
                for pc=1:size(isi_y,1)
                    [overall_dissimilarity,plot_x_values,plot_y_values(pc,:)] = SPIKY_f_pico(isi_x,isi_y(pc,:),para.tmin);
                end
                hold on
                if plot_profiles>1
                    plot(plot_x_values,plot_y_values(2:end,:))
                end
                if mod(plot_profiles,2)>0
                    plot(plot_x_values,plot_y_values(1,:),'k','LineWidth',1.5)
                end
            elseif ismember(measure,[3 4 5])                               % piecewise linear profiles (SPIKE) can be plotted right away
                x=SPIKY_loop_results.(measure_var).time;
                y=SPIKY_loop_results.(measure_var).profile;
                hold on
                if plot_profiles>1
                    plot(x,y(2:end,:))
                end
                if mod(plot_profiles,2)>0
                    plot(x,y(1,:),'k','LineWidth',1.5)
                end
            elseif measure==6                                              % PSTH
                x=SPIKY_loop_results.(measure_var).time;
                y=SPIKY_loop_results.(measure_var).profile;
                hold on
                if plot_profiles>1
                    plot(x,y(2:end,:))
                end
                if mod(plot_profiles,2)>0
                    plot(x,y(1,:),'k','LineWidth',1.5)
                end
            end
            xlim([para.tmin para.tmax])
            title ([measure_name,'   ---   Dissimilarity profile'],'FontWeight','bold','FontSize',14)
        end

        if mod(plotting,8)>3 && measure<6
            subplotc=subplotc+1;
            subplot(num_plots,1,subplotc)                                      % Dissimilarity matrix
            mat=SPIKY_loop_results.(measure_var).matrix;
            imagesc(mat)
            axis square
            if size(mat,1)<10
                set(gca,'XTick',1:size(mat,1),'YTick',1:size(mat,1))
            end
            title ([measure_name,'   ---   Dissimilarity matrix'],'FontWeight','bold','FontSize',14)
            colorbar
        end
    end
end


%%
figure;hold on;
for k=1:700
    spikes = batch{run.StimMat_str(k)};
    if ~isempty(spikes)
        scatter(spikes,ones(1,length(spikes))*k,'+');
    end
end