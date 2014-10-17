function [S,D,UPs]=findUPstates(mv, lt, ut, rmp, dt )

S = [];
D = [];
UPs = [];
combn = [];

if isempty(mv)
    return
end

% Smooth the trace first:
SR = 1000; % Sampling rate
fco = 4; % Frequency cutoff
fmv=filter(gaussfiltcoef(SR,fco),1,mv);

% fmv = [zeros(10,1);ones(10,1)*5;ones(10,1)*50;ones(10,1)*5;ones(10,1)*50;ones(10,1)*0;ones(10,1)*30;ones(10,1)*0]'-66
% Eipame me ti Yiota: Duration of depolarization above the plateau potential >>  the mean background plateau

L = (fmv-rmp)>lt ;% 2mv lower threshold , given baseline (noise) equal 1mV, as in Yousheng, 2003
U = (fmv-rmp)>ut ;% upper threshold

if size(L,2) == 1
    L = L'; % Transpose to row vector.
end
if size(U,2) == 1
    U = U'; % Transpose to row vector.
end

% locate UP states that pass upper/lower threshold
[Lm,Ls] = regexp(sprintf('%i',L),'1+','match');
[Um,Us] = regexp(sprintf('%i',U),'1+','match');

% locate the bumps that pass lower threshold and upper threshold:
for il = 1:length(Lm)
    for iu = 1:length(Um)
        if ((length(Lm{il}) >= length(Um{iu})) && (Ls(il)<=Us(iu)))
            combn(il) = true;
        end
    end
end

% keep only those:
M = {Lm{1,find(combn)}};
S = Ls(1,find(combn));
[~,I] = find(cellfun('length',M)>dt); % Duration threshold
if isempty(I)
    S = [];
    return;
else
    
    D = cellfun('length',M(I)) ; % duration of UP states
    S = S(I);  % start times of UP states
    
    % Save vm of UP states in a cell if needed:
    UPs = cell(length(D),1);
    
    for i=1:length(UPs)
        UPs{i} = mv(S(i):(S(i)-1)+D(i));
    end
    
end

end

function b=gaussfiltcoef(SR,fco)
%GAUSSFILTCOEF  Return coefficients of Gaussian lowpass filter.
% SR=sampling rate, fco=cutoff (-3dB) freq, both in Hz.
% Coeffs for FIR filter of length L (L always odd) are computed.
% This symmetric FIR filter of length L=2N+1 has delay N/SR seconds.
% Examples of use
%    Compute Gaussian filter frequency response for SR=1000, fco=50 Hz:
%    freqz(gaussfiltcoef(1000,50),1,256,1000);
%    Filter signal X sampled at 5kHz with Gaussian filter with fco=500:
%    y=filter(gaussfiltcoef(5000,500),1,X);
% SR, fco are not sanity-checked.  WCR 2006-10-11.

b=0;
a=3.011*fco;
N=ceil(0.398*SR/fco);   %filter half-width, excluding midpoint
%Width N corresponds to at least +-3 sigma which captures at least 99.75%
%of area under Normal density function. sigma=1/(a*sqrt(2pi)).
L=2*N+1;                %full length of FIR filter
for k=-N:N
    b(k+N+1)=3.011*(fco/SR)*exp(-pi*(a*k/SR)^2);
end;
%b(k) coeffs computed above will add to almost exactly unity, but not
%quite exact due to finite sampling and truncation at +- 3 sigma.
%Next line adjusts to make coeffs b(k) sum to exactly unity.
b=b/sum(b);
end