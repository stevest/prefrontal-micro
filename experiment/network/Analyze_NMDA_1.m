cd /home/stevest/Downloads/PFC_C3_5_alltogetherStefanos/experiment/network/VALIDATION_NMDA/
% cd /home/stevest/Desktop/PFC_C3_5_alltogether/experiment/network/VALIDATION_NMDA/
spikeThres=16; % was 16
eventOccure = 5000 %stimulus timming
% figure(1);
% figure(2);

Mode = '_Single';
% Mode = 'Paired';
plotColor = 'k';
for fores = 1:2
    cd(Mode);
    cnt=1;
    cnt2=1;
    
    for intensity=0:40:760
        for runs=0:0
            for currDend=0:17
                if(intensity==0)
                    intensity=1;
                end
                
                temp = sprintf('NMDA_%d_%d_%d.txt',runs,intensity,currDend)  ;
                %         temp = sprintf('NMDA_%d.txt',intensity)  ;
                if( exist(temp) )
                    
                    load(temp);
                    % check peak amplitude
                    tempPeak(cnt) = max(eval([temp(1:end-4),'(',sprintf('%d',eventOccure),':end)']) );
                    % check half-width:
                    tempEl = eval([temp(1:end-4),'(1:end)']);
                    tempEl = tempEl(5:end);% remove header info..
                    tempEl_Trimmed = tempEl(eventOccure:end);
                    [XMAX,IMAX,XMIN,IMIN] = extrema(tempEl_Trimmed);
                    %         halfAmpl = (tempEl(IMAX(2)+2700) - tempEl_Trimmed(1))/2;
                    tempZero = tempEl - tempEl_Trimmed(1);
                    
                    %             if(max(tempZero(eventOccure:end))<spikeThres)
                    halfAmpl = (max(tempEl_Trimmed) - tempEl_Trimmed(1))/2;
                    
                    if(fores==2)
                        if(tempZero(eventOccure+200)<(tempZero(eventOccure+210)))
                            START = eventOccure+200;
                        else
                            START = find(tempZero(eventOccure:end)>halfAmpl);
                            START = START(1) + eventOccure;
                        end
                    else
                        START = find(tempZero(eventOccure:end)>halfAmpl);
                        START = START(1) + eventOccure;
                    end
                    
                    
                    FINISH = find(tempZero(START+100:end)<halfAmpl);
                    FINISH = FINISH(1) + START;
                    
                    [START,FINISH]
                    
                    tempHalfWidth(cnt) = FINISH-START;
                    
                    
% %                     if isnan(tempHalfWidth(cnt))
%                         plot(tempZero,'b');hold on;
%                         plot([START,FINISH],[halfAmpl,halfAmpl],'r');hold on;
%                         pause;
% %                     end
%                     cla;
                    
                    if(isnan(tempHalfWidth(cnt)))
                        'IS NAN'
                        pause;
                        continue;
                    end
                    cnt = cnt+1;
                    
                    
                    
                    %             end %Spike thresholdnig..
                    %         [currDend intensity]
                    %         plot(tempZero,'b');hold on;
                    %             pause;
                    %             close(gcf);
                    
                end % if exists
            end
        end
        
        cnt=1;
        
        if ~isnan(mean(tempPeak)) || ~isnan(mean(tempHalfWidth))
            peakAmp(cnt2) = mean(tempPeak);
            halfWidth(cnt2) = mean(tempHalfWidth);
            tempHalfWidth=[];
            tempPeak=[];
            cnt2=cnt2+1;
        end
        
    end
%     figure(1);
%     pause
    peakAmp = peakAmp-peakAmp(1);
    peakAmp=peakAmp(1:end);
    % e=std(peakAmp,1,2);
    % errorbar(peakAmp,e,plotColor,'0');figure(gcf); hold on;
%     plot(peakAmp,plotColor);figure(gcf);hold on;
%     abscissa = linspace(1,length(halfWidth), length(halfWidth))  ;
%     scatter(abscissa,peakAmp,plotColor);figure(gcf);hold on;    
   
%     figure(2);
    halfWidth = halfWidth/10;
    halfWidth=halfWidth(1:end);
    % y=mean(halfWidth);
    % e=std(halfWidth,1,1);
    % errorbar(halfWidth,e,plotColor,'0');figure(gcf);hold on;
%     plot(halfWidth,plotColor);figure(gcf);hold on;
    
%     scatter(abscissa,halfWidth,plotColor);figure(gcf);hold on;
    

    TotalhalfWidth(fores,:) = halfWidth;
    TotalpeakAmp(fores,:) = peakAmp;
    
    Mode = '_Paired';
    plotColor = 'r';
    cd ..;
end

abscissa = linspace(1,length(TotalpeakAmp(1,1:end-2)), length(TotalpeakAmp(1,1:end-2)))  ;
scatter(abscissa,TotalpeakAmp(1,1:end-2),'k');figure(gcf);hold on;  

abscissa = linspace(1,length(TotalpeakAmp(2,1:end-2)), length(TotalpeakAmp(2,1:end-2)))  ;
scatter(abscissa,TotalpeakAmp(2,1:end-2),'r');figure(gcf);hold on;  

figure;
abscissa = linspace(1,length(TotalhalfWidth(1,1:end-2)), length(TotalhalfWidth(1,1:end-2)))  ;
scatter(abscissa,TotalhalfWidth(1,1:end-2),'k');figure(gcf);hold on;  

abscissa = linspace(1,length(TotalhalfWidth(2,1:end-2)), length(TotalhalfWidth(2,1:end-2)))  ;
scatter(abscissa,TotalhalfWidth(2,1:end-2),'r');figure(gcf);hold on;  

%%
for increment=1:18
    for intensity=10:10:60
        temp = sprintf('NMDAspikes_%d_%d',increment,intensity)  ;
        eval( ['plot(',temp,')'] );hold on;
    end
    pause;
    close all
end
