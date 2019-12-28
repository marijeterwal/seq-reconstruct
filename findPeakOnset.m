% Function to determine the onset of peaks in 1D data traces, either defined  
% as half peakheight (halfpeak method) or by a change in derivative (peakonset 
% method). First, peaks are found, then the data are smoothed. For the half 
% peakheight method the derivative of the smoothed traces is computed and a 
% buffered threshold crossing is found. 
% The function outputs the onset of the highest peak meeting the criteria as
% an index.

% The function requires two inputs: 
% - a config structure
% - 1D array with data 

% The config struct can specify the following parameters:
% * threstype: (default: 'peakonset')
% - 'halfpeak' to identify the half peakheight point and
% - 'peakonset' to determine the peak onset as defined by the parameter 'dthres'
% * peakthres: threshold to be used for peak detection (default = none)
% * nbelowthres: integer defining minimum number of values below threshold 
% required to accept threshold crossing (integer default = 3)
% * dthres: threshold for derivative change used to define the inflection
% point of the data trace, i.e. the peak onset. This value is not used for 
% halfpeak method. (double; default = tan((pi*5)/180))
% * plim: restricts the search for peaks to part of the trace. Limits are 
% specified as vector of two indices: [idx1,idx2] (array, default = none)
% * smoothwidth: width of smoothing kernel used, width in indices. For no
% smoothing, use 1 (integer, default = 5)

% The cfg defaults specified were used in the paper:
% Ter Wal et al., (under review). Human stereoEEG recordings reveal network 
% dynamics of decision-making in a rule-switching task.


% Marije ter Wal - 2019
% m.j.terwal@bham.ac.uk

function onset = findPeakOnset(cfg, data)

%% check inputs
if nargin < 2
    error('Not enough inputs, two required')
end
if sum(size(data)>1)>1 
    warning('Data contains more than 1 dimension, only the first dimension will be considered')
    data = data(:,1);
end

% check Config
if ~isfield(cfg,'threstype'); cfg.threstype = 'peakonset'; end
if ~isfield(cfg,'peakthres'); cfg.peakthres = min(data); end
if ~isfield(cfg,'nbelowthres'); cfg.nbelowthres = 3; end 
if ~isfield(cfg,'dthres'); cfg.dthres = tan((pi*5)/180); end
if ~isfield(cfg,'plim'); cfg.plim = [1,length(data)]; end
if ~isfield(cfg,'smoothwidth'); cfg.smoothwidth = 5; end


%% find the peak(s)
[~,peaks] = findpeaks(squeeze(data(cfg.plim(1):cfg.plim(2))), 'minpeakheight',cfg.peakthres);
if isempty(peaks)
    onset = NaN;
    return
end

% smooth the data
sk = 1/cfg.smoothwidth*ones(cfg.smoothwidth,1);
sdat = conv(squeeze(data), sk,'same');

breaktok = false;
peakn = 1;

% loop throught all peaks, from highest to lowest and find peaks that meet
% the criteria set by cfg. 
while ~breaktok && peakn <= length(peaks)
    %find peak within trace
    [peakheight,id] = findpeaks(squeeze(data(cfg.plim(1):cfg.plim(2))), 'npeaks',peakn, 'sortstr','descend');
    if isempty(peakheight) || length(peakheight)<peakn
        peakstart = NaN;
        peakid = NaN;
        break
    end
    peakid = id(peakn)+cfg.plim(1)-1;
    smoothclas = sdat(1:peakid);
    
    % find point left of the peak that has a derivative below threshold
    datshort = flipud(smoothclas(2:end));
    
    % two methods: 
    % 'halfpeak': find the half max point of the peak
    % 'onset': find first n values below derivative threshold
    if strcmp(cfg.threstype,'halfpeak')
        diffshort = flipud(smoothclas(2:end) < peakheight(peakn)/2);
        peakstart = find( (cumsum(diffshort)>cfg.nbelowthres & datshort<peakheight(peakn)/2), 1, 'first');
    elseif strcmp(cfg.threstype,'peakonset')
        diffshort = peakheight(peakn)*flipud(diff(smoothclas(1:end)))<cfg.dthres;
        peakstart = find( (cumsum(diffshort)>cfg.nbelowthres & diffshort>0 & datshort<peakheight(peakn)/2) ...
            | datshort<0, 1, 'first');
    end  
    
    if ~isempty(peakstart) && (peakid - peakstart) >= 1
        breaktok = true;
    else
        peakn = peakn+1;
    end
end

%% admin
if isempty(peakstart)
    peakstart = NaN;
end
if peakid - peakstart < 1
    onset = NaN;
else
    onset = peakid - peakstart;
end

end