tmp = orgW;
tmp(tmp==0) = [];
length(tmp)
ranges = 0:0.05:5;
histo_org = histcounts(tmp,ranges)/length(tmp);

tmp = newW;
tmp(tmp==0) = [];
length(tmp)
ranges = 0:0.05:5;
histo_new = histcounts(tmp,ranges)/length(tmp);

figure;hold on;
plot(ranges(1:end-1),histo_org)
plot(ranges(1:end-1),histo_new)
legend({'original','new'})
%%
tmp = orgC(1:700,1:700);
tmp = tmp .* (~eye(700));
sum(tmp(:)) / (700^2 - 700)

tmp = newC(1:700,1:700);
tmp = tmp .* (~eye(700));
sum(tmp(:)) / (700^2 - 700)



%%

tmp = W;
tmp(tmp==0) = [];

mean(tmp)
mean(W(:))

% To factor gia ta str fainetai na einai to 5:
tmp = (tmp/max(tmp(:)))*5;

ranges = 0:0.05:5;
histo = histcounts(tmp,ranges);
figure;plot(ranges(1:end-1),histo)

%%

% Giati to weight tou rnd den mou bgainei?

tmp = run.weights_str;
tmp(tmp==0) = [];

tmp = (tmp/max(tmp(:)))*2.0;

median(tmp)

ranges = 0:0.05:5;
histo = histcounts(tmp,ranges);
figure;plot(ranges(1:end-1),histo)