pathto = '/home/cluster/stefanos/Documents/GitHub/prefrontal-micro/experiment/network/DATA/'
g = dir(pathto)
files = {g(~[g.isdir]).name};
n = length(files);
rmp = -66;
distamp = [];
for fn = 1:n
    tmp = strsplit(files{fn}, {'nmda_loc_','_syn_','_nmdaw_','_ISI_','.txt'});
    id1 = str2double(tmp{2});
    id2 = str2double(tmp{3});
    id3 = str2double(tmp{4});
    id4 = str2double(tmp{5});
    vtrace = load(fullfile(pathto,files{fn}));
    distamp(id1,id2,id3,id4) = max(vtrace);
    hw(id1,id2,id3,id4) = half_width(vtrace,rmp);
end
d1 = find(~all(distamp==0,1));
d2 = find(~all(distamp==0,2));
d3 = find(~all(distamp==0,3));
d4 = find(~all(distamp==0,4));

%syn and nmda w:
da = distamp(d2,d3) - rmp;
HW = hw(d2,d3);

figure;surf(da);
set(gca,'XTickLabel',c);set(gca,'YTickLabel',r);
xlabel('synapses');ylabel('location');
axis square;

figure;surf(HW);
set(gca,'XTickLabel',c);set(gca,'YTickLabel',r);
xlabel('synapses');ylabel('location');
axis square;