function plot_ifr(st,varargin)
% PLOT_IFR(st) plots instant firing rate of spike trains given (cell
% array) in milliseconds.
% PLOT_IFR(st, tstop) plots up to tstop (in milliseconds).
% If input (st) have multiple rows, each represents one network cell.
% If input (st) have multiple columns, are ignored.
%
% author stamatiad.st@gmail.com

n = size(st,1);
if isempty(varargin)
    tstop = round(max(cell2mat(cellfun(@max,st,'uniformoutput',false))));
else
    tstop = varargin{1};
end

if n==1
    spikes = st{1,1};
    if ~isempty(spikes)
        figure;plot(instantfr(spikes,linspace(0,tstop,tstop)));
        ylabel('Instantaneous firing rate (Hz)');
        xlabel('Time (ms)');
    end
else
    IFRArray = zeros(n,tstop);
    for c=1:n
        spikes = st{c,1};
        if ~isempty(spikes)
            IFRArray(c,:) = instantfr(spikes,linspace(0,tstop,tstop));
        end
    end
    
    cm = jet(ceil(max(IFRArray(:))*100));
    cm(1,:)=[0,0,0];
    figure;imagesc(IFRArray);
    colormap(cm);
    ylabel('Cell id');
    xlabel('Time (ms)');
end



return