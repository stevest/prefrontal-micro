% blah de blah create dummy synapses:
nIncomingR=[3] % 14 (in total) synapses with the new nmda_segev
 % Background activity must result in ~1.2Hz in the post-synaptic cell as in
 % Carl C.H.Petersen Neuron, 2013
tstop = 5000
spikesRapid={};
randShiftR = 100;%180;
%rng(5,'twister');
stimPoissonR = find(poissrnd(0.1,tstop,1))  ;
    
poilR = length(stimPoissonR) ;
for c=1:1
    for i = 1:nIncomingR(1)
    spikesRapid(c,i) = { sort(round((rand(poilR,1)-0.5)*randShiftR) + stimPoissonR,'ascend')} ; %{ find(poissrnd(0.004,tstop,1)) } ; % ~12Hz
    spikesRapid{c,i}(spikesRapid{c,i}<0) = [];
    spikesRapid{c,i}(spikesRapid{c,i}>tstop) = [];
    end
end


TESTexportBackgroundStimParams(spikesRapid,mypath)