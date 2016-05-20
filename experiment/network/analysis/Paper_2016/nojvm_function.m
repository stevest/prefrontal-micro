%% Call nice function instead:
load(fullfile(osDrive(),'Documents','Glia','NetworkCreation_SN3.mat'));
%Which cluster is stimulated in each configuration:
stc_rnd = 3;
stc_str = 2;
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Rs20c%d_SN%d_spikes.mat',stc_rnd-1, run.sn)));
load(fullfile(osDrive(),'Documents','Glia','dataParsed2Matlab',sprintf('updatedStimGABAb01NEWBGST_Ss20c%d_SN%d_spikes.mat',stc_str-1,run.sn)));
run.tstop = 20000;
run.nruns = 100;

% List of stimulated/non-stimulated cells in each configuration:
sc_rnd = run.stimulatedCells_rnd{stc_rnd};
nsc_rnd = find(~ismember(1:700,run.stimulatedCells_rnd{stc_rnd}));
sc_str = run.stimulatedCells_str{stc_str};
nsc_str = find(~ismember(1:700,run.stimulatedCells_str{stc_str}));
% Array of windows Q to apply:
% Qseq = [100,50,40,30,20,10,8,6,4];
Qseq = [20,10,8,6,4];

% Na e3etazw ono to network, mono to stimulated k mono to recruited 'H ola
% ta ypolloipa kyttara (ola ane3artita apo to an einai recruited).

% Only stimulated cluster:
configuration = 'rnd';
eval( sprintf('st = batch_%s_spikes;',configuration) );
eval( sprintf('[voteState, U, smoothed_states] = createVoteState(run, Qseq, st, stc_%s-1, configuration, ''save'');',configuration) );
configuration = 'str';
eval( sprintf('st = batch_%s_spikes;',configuration) );
eval( sprintf('[~, ~, ~] = createVoteState(run, Qseq, st, stc_%s-1, configuration, ''save'');',configuration) );
