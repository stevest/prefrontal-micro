function events = nonhomogenouspoissrnd(nsamples,lambdaMax,rate)
    poissevents = poissrnd(lambdaMax,1,nsamples);
    
    nrate = rate./(max(rate)/lambdaMax);
    
%     Tevents = zeros(1,nsamples);
    thinningProb = nrate./lambdaMax;
    Tevents = logical(poissevents) & (rand(1,nsamples)<thinningProb);
    events = find(Tevents);
end

% function events = nonhomogenouspoissrnd(nsamplesM,nsamplesN,lambdaMax,rate)
%     poissevents = poissrnd(lambdaMax,nsamplesM,nsamplesN);
%     
%     nrate = rate./(max(rate)/lambdaMax);
%     
% %     Tevents = zeros(1,nsamples);
%     thinningProb = nrate./lambdaMax;
%     
%     Tevents = logical(poissevents) & (rand(nsamplesM,nsamplesN)<repmat(thinningProb,nsamplesM,1));
%     events = cell(nsamplesM,1);
%     for k=1:nsamplesM
%         events{k} = find(Tevents(k,:));
%     end
% end