function getStates(run, Qseq, st, stc, nstates, threshold, configuration)
% Returns the network states for window of size Q

% % typecheck:
% if ~isfloat(Qseq)
%     error('Error. \nInput variable Q must be a float, not a %s.',class(Qseq));
% end


% % Prepare spike arrays:
% if strcmpi(configuration,'str');
%     StimMat = run.StimMat_str;
% else
%     StimMat = run.StimMat_rnd;
% end

% labels = zeros(700,1);
% for k = 1:size(StimMat,2)
%     labels(StimMat(:,k)) = k;
% end

% NC = size(StimMat,2);
% failedToLoad = zeros(run.nPC,100,7);
% for stc=1:NC
%     for r = 1:run.nruns
%         for c = 1:run.nPC
%             failedToLoad(c,r,stc) = isempty(RUNS{1,stc}{c,r});
%         end
%     end
% end
% rurange = find(~any( squeeze(any(failedToLoad,1)),2 ));
%
% % If missing runs, remove them from the array:
% if run.nruns ~= length(rurange)
%     for stc=1:7
%         RUNS{1,stc} = RUNS{1,stc}(:,rurange);
%     end
%     run.nruns = length(rurange);
% end


m = run.nruns ; %No of chains (m)
n = run.tstop ; %No of itterations (n)
for Qs = 1:length(Qseq)
    Q = Qseq(Qs) ; % simple window (ms)
    Qr = floor(n / Q) ; % length of reshaped spiketrain array
    n_over_b = 20 ; % No of batches
    b = floor(Qr/n_over_b); % Batches' length (b) in reshaped spiketrain units
    if b < 2
        error('Batches'' length must be greater than one sample!');
    end
    
    %     for stc=1%:NC
    %         wspktrain = windowedSpikeTrain(RUNS,stc, m, n, Q, Qr);
    
    % Find unique network states:
    s = (((1:Qr)-1) * Q) + 1 ;
    e = ((1:Qr) * Q) ;
    tic;
    uniqueStates = zeros(m,Qr);
    tmparr = zeros(size(st,2),size(st,1));
    for k=1:Qr
        ss = s(k);
        ee = e(k);
        %             tic;
        % %             uniqueStates(:,k) =
        %             custombi2de(cellfun(@(x) any((ss<x)&(x<=ee)),st)');
        %             fprintf('Generating uniqueStates took (cellfun): %fs\n',toc);
        %             tic;
        
        for ii = 1:size(st,1)
            for jj = 1:size(st,2)
                tmpvar = st{ii,jj};
                tmparr(jj,ii) = any((ss<tmpvar)&(tmpvar<=ee)) ;
            end
        end
        uniqueStates(:,k) = bi2de(tmparr);
        %             fprintf('Generating uniqueStates took (loop): %fs\n',toc);
    end
    fprintf('Generating uniqueStates took: %fs\n',toc);
    
    
    %  Work with binary representations:
    %         U = unique(sort(bi2de(reshape(wspktrain,700,[])')));
    U = unique(uniqueStates(:));
    fprintf('Total No of states found: %d\n',length(U));
    % Votestate matrix and U matrix have the same ordering in states.
    voteState = zeros(size(U,1),Qr);
    % Remove stimulation period:
    for qr=ceil(1500/Q)+1:Qr
        %             B = sort(bi2de(squeeze(wspktrain(:,:,qr))'));
        [~,locStates] = ismember(uniqueStates(:,qr),U);
        voteState(:,qr) = accumarray([1:size(U,1), locStates']',[voteState(:,qr)', ones(1,length(locStates))]) ;
    end
    
%     % Sort states by maximum relative frequency:
%     [maxfreqstates,maxfreqidx] = sort(max(voteState,[],2),'descend');
%     fprintf('Most frequent State is (by max frequency) %d\n',length(U));
    %Sort by number of neurons:
%     [maxrf,idx] = sort(arrayfun(@(x) sum(de2bi(x)), U),'descend');
    % get mean relative frequency:
    voteState = voteState./m; 
    %Get most frequent states by smooth data (NOT raw relative frequency):
    tic;
    maxstates = zeros(size(voteState,1),1);
    for kk=1:size(voteState,1)
        kk
        tmp = smooth(voteState(kk,:)',100,'rlowess');
        maxstates(kk) = max(tmp(floor(Qr/2):Qr));
    end
    fprintf('Generating max of smoothed states took: %fs\n',toc);
    [maxfreqstates,maxfreqidx] = sort(maxstates,'descend') ;
    largeFreq = maxfreqstates(maxfreqstates~=0);
    idx = maxfreqidx(1:length(largeFreq));
%     % Old code (One2one 19/1/2016):
%     % Get only states with mean frequency above one threshold:
%     % Sort those states with number of active neurons:
%     largeFreq = find(mean(voteState(:,floor(Qr/2):Qr),2) > threshold) ;
%     Ul = U(largeFreq);
%     [maxrf,idx] = sort(arrayfun(@(x) sum(de2bi(x)), Ul),'descend');
%     % Get relative frequency percentage:
% %     voteState = voteState(idx,:)./m;
       
    
    
    maxy = max(maxfreqstates);
    % Old Code:
%     maxy = max(maxrf/m);

    % No of subplot columns:
    ncols=7;
    
%     idxmat = [];
%     for k=1:ceil(length(largeFreq)/(nstates))
%         idxmat = [idxmat;
%     end
    
    if true
        for s=1:ceil(length(largeFreq)/(nstates))
%             mystates = (1:nstates) + nstates*(s-1);
% ftmpind = (1:nstates) + nstates*(s-1);
            % Den to pistevw oti grafw ton parakatw kwdika:
            
            if (nstates + nstates*(s-1)) <=length(largeFreq)
                ftmpind = (1:nstates) + nstates*(s-1);
            else
                ftmpind =  floor(length(largeFreq)/nstates)*nstates+1:length(largeFreq) ;
            end
            if ncols>length(ftmpind)
                ncols = length(ftmpind);
            end
            % Get Absolute index to states array U:
            % Old code:
%             mystates = largeFreq(idx(ftmpind));
            mystates = idx(ftmpind);
            
            nrows = ceil(length(mystates)/ncols);
            figure('Name',sprintf('%s, Cluster: %d, Q: %d, From: %d TH:%f',upper(configuration),stc,Q,nstates*(s-1),threshold),'Position', [ 63           1        1858        1003]);
            for k=1:length(mystates)
                
                subplot(nrows,ncols,k);hold on;
                plot(voteState(mystates(k),:),'x');
                plot(smooth(voteState(mystates(k),:)',100,'rlowess'),'r', 'linewidth',2);
                % Set axis limits:
                xlim([0,Qr]);
                ylim([0,maxy]);
%                 title(sprintf('Active State Neurons: %d',sum(de2bi(U(idx(k)))) ));
                title(sprintf('Active State Neurons: %d',sum(de2bi(U(mystates(k)))) ));
                
                if ~(mod(k-1,ncols)==0)
                    set(gca, 'YTick', []);
                end
                if ~( (((nrows-1)*ncols)+1)<=k )
                    set(gca, 'XTick', []);
                end
                
            end
        end
    end
    
    %     end
end

end

% function d = custombi2de(b, varargin)
% %BI2DE Convert binary vectors to decimal numbers.
% %   D = BI2DE(B) converts a binary vector B to a decimal value D. When B is
% %   a matrix, the conversion is performed row-wise and the output D is a
% %   column vector of decimal values. The default orientation of the binary
% %   input is Right-MSB; the first element in B represents the least
% %   significant bit.
% %
% %   In addition to the input matrix, two optional parameters can be given:
% %
% %   D = BI2DE(...,P) converts a base P vector to a decimal value.
% %
% %   D = BI2DE(...,MSBFLAG) uses MSBFLAG to determine the input orientation.
% %   MSBFLAG has two possible values, 'right-msb' and 'left-msb'.  Giving a
% %   'right-msb' MSBFLAG does not change the function's default behavior.
% %   Giving a 'left-msb' MSBFLAG flips the input orientation such that the
% %   MSB is on the left.
% %
% %   Examples:
% %       B = [0 0 1 1; 1 0 1 0];
% %       T = [0 1 1; 2 1 0];
% %
% %       D = bi2de(B)
% %       E = bi2de(B,'left-msb')
% %       F = bi2de(T,3)
% %
% %   See also DE2BI.
%
% %   Copyright 1996-2011 The MathWorks, Inc.
% % SSS removed checks rof speed
%
% b = double(b);  % To allow non-doubles to work
%
% % Set up the base to convert from.
% p = 2;
%
% n = size(b,2);
%
% %%% The conversion
% max_length = 1024;
% pow2vector = p.^(0:1:(n-1));
% size_B = min(max_length,n);
% d = b(:,1:size_B)*pow2vector(:,1:size_B).';
%
% end

function wspktrain = windowedSpikeTrain(st,stc, m, n, Q, Qr)
% Calculate spike train in time windows

s = (((1:Qr)-1) * Q) + 1 ;
e = ((1:Qr) * Q) ;

% Number of neurons:
nlen = size(st{1,stc},1);
wspktrain = zeros(nlen,m,Qr);

for ru = 1:m
    spktrain = zeros(nlen,n);
    for c=1:nlen % must be row vector!!
        if ~isempty(st{1,stc}{c,ru})
            spktrain(c,round(st{1,stc}{c,ru}.spikes')) = 1;
        end
    end
    for k=1:Qr
        wspktrain(:,ru,k) = any(spktrain(:,s(k):e(k)),2) ;
    end
end
end