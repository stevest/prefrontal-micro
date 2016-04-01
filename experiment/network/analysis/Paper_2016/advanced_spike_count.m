
function  [number_of_spikes, spike_timing]= advanced_spike_count(vt, lt, ut, varargin)
% ADVANCED_SPIKE_COUNT return spike timings of voltage trace.
%   ADVANCED_SPIKE_COUNT(vt, lt, ht) find spikes in voltage trace vt ( that
%   first cross high threshold and again low threshold).
%   ADVANCED_SPIKE_COUNT(vt, lt, ht, 'plot') also plot results.
%
%   This updated function can handle:
%   > voltage train without spikes
%   > noisy (high freq) voltage train
%
%   author stamatiad.st@gmail.com

p = false;
for a=1:length(varargin)
    if strcmp(varargin{a},'plot');p=true;end
end

if isempty(vt)
    %You can not return number of spikes of an empty vt!
    warning('Empty voltage trace as input!');
    number_of_spikes=[];
    spike_timing=[];
    return;
end
spikeRegions = zeros(1,length(vt));

%Find values above high threshold:
a=find(vt>ut);
for k=1:length(a)
    vpart=a(k);
    if vpart==1
        continue
    elseif vt(vpart-1)<ut
        %Find were spike region ends (old code might produce
        %multiple spike regions, each producing a max and essentially
        %producing extra spikes).
        b=find(vt((vpart+1):length(vt))<lt, 1, 'first');
        % Register spike region:
        spikeRegions(vpart:vpart+b) = 1;
    end
end

[SpkRegs,spikeRegionStartIdx] = regexp(sprintf('%i',spikeRegions),'1+','match');

if isempty(SpkRegs)
    number_of_spikes=0;
    spike_timing=[];
else
    % Locate spikes by getting the max of each registered spike region:
    spikeRegionEndIdx = spikeRegionStartIdx + cellfun(@length, SpkRegs)-1;
    voltageSpikesRegions = cellfun(@(s,e) vt(s:e),...
        mat2cell(spikeRegionStartIdx,[1],ones(1,length(spikeRegionStartIdx))),...
        mat2cell(spikeRegionEndIdx,[1],ones(1,length(spikeRegionStartIdx))),'uniformoutput',false );
    [~,loc] = cellfun(@(x) max(x), voltageSpikesRegions);
    spike_timing=(spikeRegionStartIdx + loc-1);
    number_of_spikes=size(spike_timing,2);
end

% Plot results:
if (p)
    figure;hold on;
    plot(vt);
    vmax=max(vt);vmax=vmax+vmax*0.1;
    scatter(spike_timing,ones(1,number_of_spikes)*vmax,'rv','fill');
    hlt=plot([0,length(vt)],[lt,lt],'--g');
    hht=plot([0,length(vt)],[ut,ut],'--c');
    legend([hlt,hht],{'Low Threshold','High Threshold'},'Location','southeast')
end

return
