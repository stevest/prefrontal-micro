function [st,tstop]=load_raw_batch(pathto,varargin)
% Load experiment batch..
% ...GIVEN that one run is inside the pathto folder.
% One run is equivalent to one MPI batch job (in its own folder!).


s=false;
for a=1:length(varargin)
    if strcmp(varargin{a},'spikes');s=true;end
end

% jungle
r = dir(pathto);
rf = {r(~[r.isdir]).name};
% Get only .bin files:
files = sort(rf(cellfun(@(x) strcmp(x(end-3:end),'.bin'),rf)));
n = size(files,2);

if (isempty(files))
   warning('No files found in folder!');
   st={};
   return;
end
st = cell(n,1);

textprogressbar('Loading batch: ');
for fn=1:n
    %These, tmp, index change according to dataset!
    tmp = strsplit(files{fn}, {'_','.'});
    % get cell id:
    cid = str2double(tmp{1})+1;
    % get run id
    rid = str2double(tmp{2})+1;
    
    filename = fullfile(pathto,files{fn});
%     disp(filename);
    if( exist(filename,'file') )
        try
            trace = nrn_vread(filename,'n');
        catch e
            warning('Error with nrn_vread() !');
        end
        % dt equals 0.1:
        trace = trace(1:10:end);
        if s
            [~, st{cid,1}] = advanced_spike_count(trace,-20,0);
        else
            st{cid,1} = trace';
        end
    else
        st{cid,1} = [];
    end
    textprogressbar((fn/n)*100);
end
textprogressbar('done');

%given that batch will have the same tstop:
tstop = size(trace,1)-1;

return

function textprogressbar(c)
% This function creates a text progress bar. It should be called with a 
% STRING argument to initialize and terminate. Otherwise the number correspoding 
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

% Main 

if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];  
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
