load('soma_run0.dat');
[B,IDX] =sort(soma_run0(:,1));
SORTED = soma_run0(IDX,:);

for i=1:30
    figure;plot(SORTED(i,:))
end



%%
close all; clear all; clc;  