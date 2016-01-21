load(['multi_core/bg_dend_0.txt']);
load(['multi_core/bg_dend_1.txt']);
load(['multi_core/bg_dend_2.txt']);
load(['multi_core/bg_dend_3.txt']);
load(['multi_core/bg_dend_4.txt']);

figure();hold on;
scatter(bg_dend_0,ones(1,32)*1,30,rand(1,3))
scatter(bg_dend_1,ones(1,32)*2,30,rand(1,3))
scatter(bg_dend_2,ones(1,32)*3,30,rand(1,3))
scatter(bg_dend_3,ones(1,32)*4,30,rand(1,3))
scatter(bg_dend_4,ones(1,32)*5,30,rand(1,3))