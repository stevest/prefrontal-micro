function [S,D,UPs]=findUPstates(mv, lt, ut, rmp, dt )
    
% No convolution implemented yet...
%     % Smooth the trace first:
%     mv = RUNS_str{1,t}{c,ru}.mv;
%     mv = mv - mean(mv) ; 
%     plot(mv);hold on;
%     gf = gausswin(500) ;
%     gf = gf/norm(gf) ;
%     mv = conv(mv, gf,'same') ;
%     plot(mv,'r');
    
    A = (mv-rmp)>lt ;% 2mv lower threshold , given baseline (noise) equal 1mV, as in Yousheng, 2003
    % upper threshold not implemented yet...

    if size(A,2) == 1
        A = A'; % Transpose to row vector.
    end

    % locate UP states that pass lower threshold
    [M,S] = regexp(sprintf('%i',A),'1+','match'); 
    [~,I] = find(cellfun('length',M)>dt) % Duration threshold
    D = cellfun('length',M(I)) ; % duration of UP states
    S = S(I);  % start times of UP states
    
    % Save vm of UP states in a cell if needed:
    UPs = cell(length(D),1);
    
    for i=1:length(UPs)
        UPs{i} = mv(S(i):S(i)+D(i)-1);
    end

end