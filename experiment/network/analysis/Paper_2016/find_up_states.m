function [S,D,UPs]=find_up_states(vt, lt, ut, rmp, durt, varargin )
% FIND_UP_STATES detects UP states in voltage trace.
% FIND_UP_STATES(vt, lt, ht, rmp, durt) returns UP states in the vt(mV)
% voltage trace defined by lower and upper thresholds lt(mV) ut(mV) as in:
% Shu et. all,Barrages of synaptic activity control the gain and
% sensitivity of cortical neurons., 2003
%
% Duration threshold durt(ms) can be applied to filter fast UP states.
%
% Returns the Starting index S of each UP state in vt and its duration D.
% Voltage trace of each UP state is contained in UPs.
%
%   author stamatiad.st@gmail.com

S = [];
D = [];
UPs = {};
ah = nan;
p = false;
for a=1:length(varargin)
    if strcmp(varargin{a},'plot');p=true;end
    if isa(varargin{a},'matlab.graphics.axis.Axes');ah = varargin{a};end
end


if isempty(vt)
    return
end


% Smooth the trace first:
% SR = 1000; % Sampling rate
% fco = 4; % Frequency cutoff
% fvt=filter(gaussfiltcoef(SR,fco),1,vt);
fvt=smooth(vt,100,'moving')';
if 1%ampt < ut
    L = (fvt-rmp)>lt ;% 2mv lower threshold , given baseline (noise) equal 1mV, as in Yousheng, 2003
    U = (fvt-rmp)>ut ;% upper threshold
    
    % locate UP states that pass upper/lower threshold
    M = sum([U;L]);
    [Ml,Ms] = regexp(sprintf('%i',L),'1+','match');
    Md = cellfun('length',Ml);
    
    I = false(1,length(Ml));
    for k=1:length(Ms)
        I(k) = any(M(Ms(k):Ms(k)+Md(k)-1)==2) && (Md(k)>=durt);
    end
    
    % Naively remove the first UP state (smooth artifact):
    % if fvt(1)>lt+rmp;I(1) = false;end
    
    if isempty(I)
        S = [];
        return;
    else
        D = Md(I);
        S = Ms(I);  % start times of UP states
        % Save vm of UP states in a cell if needed:
        UPs = cell(length(D),1);
        for i=1:length(UPs)
            UPs{i} = vt(S(i):(S(i)-1)+D(i));
        end
    end
    
end

if(p)
    fUPs = cell(length(D),1);
    for i=1:length(fUPs)
        fUPs{i} = fvt(S(i):(S(i)-1)+D(i));
    end
    if ~isa(ah,'matlab.graphics.axis.Axes')
        figure;
        ah = axes;hold on;
    end
    
    plot(ah,vt,'color',[176/256,195/256,247/256])
    for k = 1:length(D)
        ha = area(ah,S(k):S(k)+D(k)-1,fUPs{k},rmp,'LineStyle','none');
        ha.FaceColor = [34/256,115/256,237/256];
%         ha.FaceAlpha = 0.5;
    end
    hlt = plot(ah,[0,length(vt)],[lt+rmp,lt+rmp],'--','color',[237/256,154/256,34/256]);
    hut = plot(ah,[0,length(vt)],[ut+rmp,ut+rmp],'--','color',[237/256,68/256,34/256]);
    
    plot(fvt,'color',[19/256,133/256,54/256])
    legend([hlt,hut], {sprintf('Low Threshold: %.3f',lt),sprintf('Upper Threshold: %.3f',ut)});
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