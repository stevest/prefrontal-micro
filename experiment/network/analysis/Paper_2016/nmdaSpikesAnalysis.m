pathto = '/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/DATA/'
g = dir(pathto);
files = {g(~[g.isdir]).name};
n = length(files);
rmp = -66;
distamp = [];
vtrace={};
hw = [];
for fn = 1:n
    tmp = strsplit(files{fn}, {'nmda_loc_','_syn_','_nmdaw_','_ISI_','.txt'});
    id1 = str2double(tmp{2});
    id2 = str2double(tmp{3});
    id3 = str2double(tmp{4});
    id4 = str2double(tmp{5});
    vtrace{id1,id2,id3,id4} = load(fullfile(pathto,files{fn}));
    distamp(id1,id2,id3,id4) = max(vtrace{id1,id2,id3,id4});
    hw(id1,id2,id3,id4) = half_width(vtrace{id1,id2,id3,id4},rmp);
end

save('NMDA_Bistability.mat','vtrace','distamp','hw','-v7.3');

% % get id for non empty elements:
% % not working:
% d1 = find(~all(distamp==0,1));
% d2 = find(~all(distamp==0,2));
% d3 = find(~all(distamp==0,3));
% d4 = find(~all(distamp==0,4));

distamp2 = squeeze(distamp(70,:,:,10));
hw2 = squeeze(hw(70,:,:,10));
d1 = find(~all(distamp2==0,1));
d2 = find(~all(distamp2==0,2));
distamp3 = distamp2(d2,d1);
hw3 = hw2(d2,d1);

%syn and nmda w:
max_amplitude = distamp3 - rmp;
HW = hw3;

figure;surf(max_amplitude,'linestyle','none');
title('Max Amplitude');
set(gca,'XTickLabel',d1);set(gca,'YTickLabel',d2);
xlabel('synapses');ylabel('nmda weight');
axis square;

figure;surf(HW,'linestyle','none');
title('Half width');
set(gca,'XTickLabel',d1);set(gca,'YTickLabel',d2);
xlabel('synapses');ylabel('nmda weight');
axis square;