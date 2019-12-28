% Multidimensional scaling to project a 2D distance matrix onto a single axis. 
% This code generates a 1D axis and iteratively updates the locations of
% the N items on this axis to reduce the difference with the distance matrix.
% The code stops when the order of the items does not change anymore or
% when a maximum number of iterations is reached.

% The function takes the following inputs: 
% - a config struct that contains the parameters (see below)
% - the N x N distance matrix Mdist, where N is the number of items in the dataset
% - an optional weight matrix Mweight, also N x N
% Usage without weight matrix: 
% [timesOut,locOut,updateOut] = mds1d(cfg,Mdist,[])

% The cfg struct can contain the following fields:
% * maxiter: the maximum number of iterations performed. If the code
% converges before reaching maxiter, it will stop before this number
% (default = 1000);
% * learnpar: multiplication factor that determines the weight of the updates
% on consecutive iterations (default = 0.9);
% * updateThres: integer that specifies the number of consecutive
% iterations without change that is required to accept the result. When
% this number of iterations is reached, the code stops updating the
% x-values. To always use maxiter iterations, set cfg.updateThres = 0
% (default = 10);
% * plotOnline: switch to plot the iterations while they are computed (default = false);
% * plot: switch to plot the end result of the MSD procedure (default =
% false).

% The function produces the following outputs:
% - timesOut: vector with the x-values for the N items
% - locOut: vector with the final order of the N items
% - updateOut: Ni x N matrix containing all intermediate orderings of the N
% items for each of the Ni iteration

% To see the code in action, run:
%{
cfg = []; cfg.plotOnline = true;
Dist = rand(10);
[timesOut,locOut,updateOut] = mds1d(cfg,0.2*Dist-0.1);
%}

function [timesOut,locOut,updateOut] = mds1d(cfg,Mdist,Mweight)

%% check config & inputs

if nargin < 2
    error('Not enough input arguments')
end

if nargin == 2
    Mweight = [];
end

if ~isfield(cfg,'maxiter'); cfg.maxiter = 1e3; end
if ~isfield(cfg,'learnpar'); cfg.learnpar = 0.9; end
if ~isfield(cfg,'updateThres'); cfg.updateThres = 10; end 
if ~isfield(cfg,'plotOnline'); cfg.plotOnline = false; end
if ~isfield(cfg,'plot'); cfg.plot = false; end

if isempty(Mweight)
    Mweight = ones(size(Mdist));
end

%% initialize 

% define an (arbitrary) initial time line
nloc = size(Mdist,1);
tspan = linspace(min(Mdist(:)),max(Mdist(:)),2*nloc);%linspace(-0.5,0.15,66);

% initial state is random
regTime = tspan(randperm(length(tspan),nloc));

if cfg.plotOnline
    reglabels = cellstr(num2str([1:nloc]'));
    
    f1 = figure; hold on;
    plot(regTime(1,:),ones(nloc,1),'o-')
    text(regTime(1,:),1.1*ones(nloc,1), reglabels,'fontsize',6)
    xlabel('Relative time (s)')
    set(gca,'ytick',[]);
    f2 = figure; hold on;
elseif cfg.plot
    f2 = figure; hold on;
end

%% update locations 

[~,updateOut(1,:)] = sort(regTime);
updateCheck = 0;
ni = 2;

while updateCheck<=cfg.updateThres && ni<=cfg.maxiter
    
    % update all but one locations
    for r = 2:nloc
        
        % get neighbours
        nb = setdiff(1:nloc,r); 
        dt = zeros(length(nb),1);
        
        for rr = 1:length(nb)
            % weight by variance
            dd = Mweight(r,nb(rr));
            % error in time difference = 
            % [current time difference] - [desired time difference]
            d = (regTime(r) - regTime(nb(rr))) - Mdist(r,nb(rr));
            
            % find updates based on each neighbour
            if ~isnan(d)
                dt(rr) = sign(d)*cfg.learnpar*(exp(abs(d*dd))-1);
            else
                dt(rr) = NaN;
            end
        end
        
        % do average update for location r
        update = nanmean(dt);
        regTime(r) = regTime(r) - update;
    end
    [~,updateOut(ni,:)] = sort(regTime);
    
    % plot online
    if cfg.plotOnline
        figure(f1); clf
        plot(regTime,ones(nloc,1),'o-')
        text(regTime,1.1*ones(nloc,1), reglabels,'fontsize',6)
        xlabel('Relative time (s)')
        set(gca,'ytick',[]);
        pause(0.1)
    end
    
    % update the loop variables
    if nansum(abs(updateOut(ni,:)-updateOut(ni-1,:))) == 0
        updateCheck = updateCheck + 1;
    else
        updateCheck = 0;
    end
    ni = ni+1;
end

%% admin
[sorttimes,sID] = sort(regTime);
timesOut = sorttimes-sorttimes(1);
locOut = sID;

if cfg.plotOnline || cfg.plot
    reglabels = cellstr(num2str([1:nloc]'));
    figure(f2);
    clf
    hold on
    plot(timesOut,1:nloc,'o-')
    text(timesOut-0.005,[1:nloc]+0.5, reglabels(sID),'fontsize',6)
    xlim([-0.02,max(timesOut)-timesOut(1)+0.02])
    xlabel('Relative time (s)')
    ylabel('Order')
end

end