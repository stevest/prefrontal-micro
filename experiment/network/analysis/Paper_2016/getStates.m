function [voteState, U] = getStates(run, Q, st)
% Returns the network states for window of size Q

% % typecheck:
% if ~isfloat(Qseq)
%     error('Error. \nInput variable Q must be a float, not a %s.',class(Qseq));
% end

% % If missing runs, remove them from the array:
% if run.nruns ~= length(rurange)
%     for stc=1:7
%         RUNS{1,stc} = RUNS{1,stc}(:,rurange);
%     end
%     run.nruns = length(rurange);
% end

fprintf('Initializing getStates (heavy)...\n');
m = size(st,2) ; %No of chains (m)
n = run.tstop ; %No of itterations (n)
N = size(st,1);
% Q is the simple window (ms)
Qr = floor(n / Q) ; % length of reshaped spiketrain array
n_over_b = 20 ; % No of batches
b = floor(Qr/n_over_b); % Batches' length (b) in reshaped spiketrain units
if b == 1
    warning('Batch is equal to Q!');
end
if b < 1
    error('Batches'' length must be greater than one sample!');
end

fprintf('Generating wspktrain...');
% wspktrain = zeros(N,m,Qr);
tic;

% sleek parfor
wst = reshape(st,1,[])';
winst = cell(m*N,1);
parfor c = 1:m*N
    if ~isempty(wst{c})
        tmpst = zeros(1,n);
        tmpwinst = zeros(1,Qr);
        tmpst(round(wst{c})) = 1;
        for k=1:Qr
            tmpwinst(k) = any( tmpst( ((((k)-1)*Q)+1):((k) * Q) ) ,2) ;
        end
       winst{c} =  tmpwinst;
    else
        winst{c} = zeros(1,Qr);
    end
end

% wstates = reshape(winst,N,[])';
wspktrain =  cell2mat( cellfun(@(x) reshape(x,1,1,[]),reshape(winst,N,[]),'uniformoutput',false) );
clear wst winst;

fprintf(' took: %fs\n',toc);
%Work with string representations:
% allStates = mat2cell( reshape(wspktrain,N,[])', ones(m*Qr,1),N );
fprintf('Getting unique states...');
acrossQs = squeeze(mat2cell( permute(wspktrain,[2,1,3]), m,N,ones(Qr,1) ));
strAcrossQs = cell(length(acrossQs),1);
for k=1:length(acrossQs)
    for j = 1:m
        strAcrossQs{k,1}{j,1} = sprintf('%i',acrossQs{k}(j,:));
    end
end

A = cellfun(@(x) unique(x), strAcrossQs,'uniformoutput',false);
U = unique(vertcat( A{:} ));
clear A;

fprintf('took: %fs\n',toc);
% voteState = zeros(size(U,1),Qr);
voteState = cell(Qr,1);
for k=1:Qr
    voteState{k} = zeros(length(U),1);
end
fprintf('Generating voteStates...');
tic;
parfor qr=1:Qr
    [~,locStates] = ismember(strAcrossQs{qr},U);
%     if any(locStates==0)
%         error('This should never happen.. Location can not be zero.');
%     end
    voteState{qr} = accumarray([1:size(U,1), locStates']',[voteState{qr}', ones(1,length(locStates))]) ;
end
voteState = cell2mat(voteState');
fprintf('took: %fs\n',toc);

% validate vote state matrix:
if ~all(sum(voteState) == m)
    error('voteState validation failed!');
end
fprintf('Terminating getStates... (phew!)\n');

return;
