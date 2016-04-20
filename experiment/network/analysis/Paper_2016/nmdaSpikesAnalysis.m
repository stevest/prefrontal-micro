pathto = 'C:\Users\user\Documents\DATA2';
g = dir(pathto);
files = {g(~[g.isdir]).name};
n = length(files);
rmp = -66;
distamp = [];
vtrace={};
hw = [];
for fn = 1:n
    fn
    tmp = strsplit(files{fn}, {'nmda_loc_','_syn_','_nmdaw_','_ISI_','.txt'});
    id1 = str2double(tmp{2});
    id2 = str2double(tmp{3});
    id3 = str2double(tmp{4});
    id4 = str2double(tmp{5});
    vtrace{id2,id3} = load(fullfile(pathto,files{fn}));
    distamp(id2,id3) = max(vtrace{id2,id3});
    hw(id2,id3) = half_width(vtrace{id2,id3},rmp);
end

% save('C:\Users\user\Documents\NMDA_Bistability_v6.mat','vtrace','distamp','hw','-v6');

% % get id for non empty elements:
% % not working:
% d1 = find(~all(distamp==0,1));
% d2 = find(~all(distamp==0,2));
% d3 = find(~all(distamp==0,3));
% d4 = find(~all(distamp==0,4));

% distamp2 = squeeze(distamp(70,:,:,10));
% hw2 = squeeze(hw(70,:,:,10));
% d1 = find(~all(distamp2==0,1));
% d2 = find(~all(distamp2==0,2));
% distamp3 = distamp2(d2,d1);
% hw3 = hw2(d2,d1);

d1 = find(~all(distamp==0,2));
d2 = find(~all(distamp==0,1));

%syn and nmda w:
max_amplitude = distamp(d1,d2) - rmp;
max_amplitude = max_amplitude(1:end-1,:)' ;
HW = hw(d1,d2);
VM = vtrace(d1,d2);

x = 1:length(d1)-1;
y = 1:length(d2);
synsTick = 1:5:length(d1);
nmdawTick = 1:10:length(d2);
syns = d1(synsTick);
nmdaw = d2(nmdawTick)/1000;
figure;surf(max_amplitude,'linestyle','none');
title('Max Amplitude');
set(gca,'XTick',synsTick);set(gca,'YTick',nmdawTick);
set(gca,'XTickLabel',syns);set(gca,'YTickLabel',nmdaw);
xlabel('synapses');ylabel('nmda weight');
axis square;

figure;surf(HW,'linestyle','none');
title('Half width');
set(gca,'XTickLabel',d1);set(gca,'YTickLabel',d2);
xlabel('synapses');ylabel('nmda weight');
axis square;

% analyze shape:
