% clear all; close all;clc;
format long g
for i=1:23
        DIR = sprintf('/home/stevest/Downloads/PFC_C3_5_alltogetherStefanos/experiment/network/VALIDATION_NMDA_K/%d/',1) ;
        cd(DIR);
        eventOccur = 5000;
        % validImg = imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/NMDA.jpg']);
        % figurel;imagesc(validImg);hold on;
        
        dataSTR = sprintf('NMDA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
        hold on;
        accumRun(i,1) = eval( ['min(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
end
mean(accumRun)


%% NMDA
for i=1:17
        DIR = sprintf('/home/stevest/Downloads/PFC_C3_5_alltogetherStefanos/experiment/network/VALIDATION_NMDA_K/%d/',1) ;
        cd(DIR);
        eventOccur = 5000;
        % validImg = imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/NMDA.jpg']);
        % figurel;imagesc(validImg);hold on;
        
        dataSTR = sprintf('NMDA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
        hold on;
        accumRun(i,1) = eval( ['max(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
end

for i=1:40
        DIR = sprintf('/home/stevest/Downloads/PFC_C3_5_alltogether/experiment/network/VALIDATION_NMDA_K/%d/',1) ;
        cd(DIR);
        eventOccur = 5000;
        % validImg = imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/NMDA.jpg']);
        % figurel;imagesc(validImg);hold on;
        
        dataSTR = sprintf('NMDA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
        hold on;
        accumRun(i,2) = eval( ['max(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
    
end
mean(accumRun')

%% NMDA Kinetics:
runsTotal = 51   ;
% runsTotal = 240   ;
cc=hsv(runsTotal);
for i=1:runsTotal
        DIR = sprintf('/home/stevest/Downloads/PFC_C3_5_alltogetherStefanos/experiment/network/VALIDATION_NMDA_K/%d/',1) ;
        cd(DIR);
        eventOccur = 5000;
        dataSTR = sprintf('NMDA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ',''Color'',',sprintf('cc(%d,:)',i), ')' ] );
        eval( ['KINETICS(',sprintf('%d',i),',:)=' dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),');'] );
        hold on;
        accumRun(i,1) = eval( ['max(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
end
mean(accumRun)

 % Plot over recordings:
 validImg = rgb2gray(imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/NMDA.jpg']));
 validImg = validImg(end:-1:1,:);
 
 imgOrdn = size(validImg,1);
 imgAbsc = size(validImg,2);
 scaleOrdn = 0.06;
 scaleAbsc = 3000;
 
 scaleFactorY = imgOrdn / scaleOrdn;
 
 KINETICS  = KINETICS * scaleFactorY;
 
 scaleFactorX = imgAbsc / scaleAbsc;
 
 Abscissa = 0:scaleFactorX:(scaleAbsc*scaleFactorX);
 
 figure;imagesc(validImg);hold on;colormap(gray);
 set(gca,'YDir','normal');
 for i=1:runsTotal
 handleMe = plot(Abscissa, KINETICS(i,:),'Color', cc(i,:), 'lineWidth', 3); hold on;
%  pause;
%  delete(handleMe);
 end

 
%% AMPA:
for i=1:1
    for run=2
        DIR = sprintf('/home/stevest/Downloads/PFC_C3_5_alltogether/experiment/network/VALIDATION_AMPA/%d/',run) ;
        cd(DIR);
        eventOccur = 5000;
        % validImg = imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/AMPA.jpg']);
        % figure;imagesc(validImg);hold on;
        
        dataSTR = sprintf('AMPA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+-2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
        hold on;
        accumRun(i,run) = eval( ['min(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
    end
    
end


% mean(accumRun')
%% AMPA Kinetics:
runsTotal = 50   ;
cc=hsv(runsTotal);
for i=1:runsTotal
        DIR = sprintf('/home/stevest/Desktop/PFC_C3_5_alltogether/experiment/network/VALIDATION_AMPA_K/%d/',1) ;
        cd(DIR);
        eventOccur = 5000;
        dataSTR = sprintf('AMPA_%d.txt',i);
        load(dataSTR);
        eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
        eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ',''Color'',',sprintf('cc(%d,:)',i), ')' ] );
        eval( ['KINETICS(',sprintf('%d',i),',:)=' dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),');'] );
        hold on;
        accumRun(i,1) = eval( ['min(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
end
mean(accumRun)


 % Plot over recordings:
 validImg = rgb2gray(imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/AMPA.jpg']));
 
 imgOrdn = size(validImg,1);
 imgAbsc = size(validImg,2);
 scaleOrdn = 0.06;
 scaleAbsc = 3000;
 
 scaleFactorY = imgOrdn / scaleOrdn;
 
 KINETICS  = KINETICS * scaleFactorY;
 
 scaleFactorX = imgAbsc / scaleAbsc;
 
 Abscissa = 0:scaleFactorX:(scaleAbsc*scaleFactorX);
 
 figure; imagesc([0,imgAbsc],[0,-imgOrdn],validImg);hold on;colormap(gray);
 set(gca,'YDir','normal');
 for i=1:runsTotal
 handleMe = plot(Abscissa, KINETICS(i,:),'Color', cc(i,:), 'lineWidth', 3); hold on;
%  pause;
%  delete(handleMe);
 end

%%
mean(accumRun')


for run=1:10
DIR = sprintf('/home/stevest/Desktop/PFC_C3_5_alltogether/experiment/network/VALIDATION_NMDA/%d/',run) ;
cd(DIR);
eventOccur = 5000;
% validImg = imread(['/home/stevest/Dropbox/Public/Neuroscience/Poirazi/NMDA.jpg']);
% imagesc(validImg);hold on;
for i=1:20
    dataSTR = sprintf('NMDA_%d.txt',i);
    load(dataSTR);
    eval( [dataSTR(1:end-4),'=',dataSTR(1:end-4),'-',dataSTR(1:end-4),'(',sprintf('%d',eventOccur+2),');' ] );
%     eval( ['plot(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
%     axis([0 3000 0 0.06 ]);hold on;
%     pause;
%     cla(gca);
    accumRun(i,run) = eval( ['max(', dataSTR(1:end-4),'(',sprintf('%d:%d',eventOccur,eventOccur+3000),')' , ')' ] );
end

end


mean(accumRun')



%%

figure('Color', [1,1,1]);
% img = imread('AMPA.jpg');
% imagesc(img);hold on;
plot(AMPA_ALPHA_100x2E00_BETA_00x2E61(490/0.1:540/0.1),'g');hold on;
axis([ 0 500 -0.1 -0.04 ])
set(gcf,'Position',[102  31  1024   666]);

%% NMDA
plot(NMDA_ALPHA_100x2E00_BETA_00x2E61(200:end),'g');hold on;
%%
figure('Color', [1,1,1]);
% img = imread('AMPA.jpg');
% imagesc(img);hold on;
plot(NMDA_ALPHA_100x2E00_BETA_00x2E61(290/0.1:1340/0.1),'g');hold on;
axis([ 0 500 -0.1 -0.04 ])
set(gcf,'Position',[102  31  1024   666]);
%%
plot(NMDA_TAU2_100(200:end),'r');hold on;
plot(NMDA_TAU2_110(200:end));hold on;
plot(NMDA_TAU2_120(200:end));hold on;
plot(NMDA_TAU2_130(200:end));hold on;
plot(NMDA_TAU2_140(200:end));hold on;
plot(NMDA_TAU2_150(200:end));hold on;
plot(NMDA_TAU2_160(200:end));hold on;
plot(NMDA_TAU2_170(200:end));hold on;
plot(NMDA_TAU2_180(200:end));hold on;
plot(NMDA_TAU2_190(200:end));hold on;
plot(NMDA_TAU2_200(200:end));hold on;
%% Tau1 variable
plot(NMDA_TAU1_50x2E000(3950/0.1:4200/0.1),'r');hold on;
plot(NMDA_TAU1_50x2E050(3950/0.1:4200/0.1),'g');hold on;
plot(NMDA_TAU1_50x2E100(3950/0.1:4200/0.1),'b');hold on;
plot(NMDA_TAU1_50x2E150(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_50x2E200(3950/0.1:4200/0.1),'y');hold on;
plot(NMDA_TAU1_50x2E250(3950/0.1:4200/0.1),'m');hold on;
%% Tau1 static
for i=1:14
load(['NMDA_TAU1_',sprintf('%d',i),'.txt']);
end


plot(NMDA_TAU1_1(3950/0.1:4200/0.1),'r');hold on;
plot(NMDA_TAU1_2(3950/0.1:4200/0.1),'g');hold on;
plot(NMDA_TAU1_3(3950/0.1:4200/0.1),'b');hold on;
plot(NMDA_TAU1_4(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_5(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_6(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_7(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_8(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_9(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_10(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_11(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_12(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_13(3950/0.1:4200/0.1),'c');hold on;
plot(NMDA_TAU1_14(3950/0.1:4200/0.1),'c');hold on;

% extract mean
MEAN = [...
max(NMDA_TAU1_1(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_2(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_3(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_4(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_5(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_6(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_7(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_8(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_9(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_10(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_11(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_12(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_13(3950/0.1:4200/0.1)), ...
max(NMDA_TAU1_14(3950/0.1:4200/0.1)) ] ;
mean(MEAN)-96.045

%%
for i=1:14
load(['AMPA_',sprintf('%d',i),'.txt']);
end

for i=1:14
eval(  ['plot(AMPA_',sprintf('%d',i),'(490/0.1:540/0.1));hold on;'] );
end

% extract mean
MAXX = [];
for i=1:14
    eval(  ['MAXX(i) =  min(AMPA_',sprintf('%d',i),'(490/0.1:540/0.1));'] );
end

mean(MAXX)- (-0.04952)




