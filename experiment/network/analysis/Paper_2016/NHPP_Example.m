snormal = @(x,m,a,sigma) ((1+erf((a.*x)/sqrt(2)))/2) .* normpdf(x,m,sigma);
kernel = snormal(-10:0.01:10,0,10,0.9);
tstop = 2000;
nevents = 3;

rate = zeros(1,tstop);

% afto mporei na proerxetai k apo poisson, opws eipes.
eventLocations = [300, 400, 700];

preconv = zeros(1,tstop);
preconv(eventLocations) = 10;
    
rate = conv(preconv(1,1:tstop),kernel,'same');
% tevents = nonhomogenouspoissrnd(2000,0.12,rate);
nsamples = 2000;
lambdaMax = 0.12;

poissevents = poissrnd(lambdaMax,1,nsamples);
poissevents(poissevents>1) = 1;

nrate = rate./(max(rate)/lambdaMax);
% Accept event at t, with probability ? / ?max
thinningProb = nrate./lambdaMax;
Tevents = logical(poissevents) & (rand(1,nsamples)<thinningProb);
events = find(Tevents);

figure;hold on;
plot([0,tstop],[lambdaMax,lambdaMax],'--');
plot(nrate)
plot(thinningProb)
scatter(find(poissevents), zeros(1,length(find(poissevents))),'k+')
scatter(events, zeros(1,length(events)),'or')
legend({'?max','?','Thinning P','All Events (?max)','Thinned Events (?(t))'})
