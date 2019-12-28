% Assessment of clustering of rankings by comparing the euclidean distance
% of the ranks in the dataset with those of a random rank dataset.

% Inputs:
% - ranks: a M x N matrix of ranks, were M is the number of runs (bootstraps,
% methods, ...) and N is the number of regions of interest
% - nperm: the number of permutations for the random dataset

% Outputs: 
% - dZ: Z-score of the rank distribution of the data against the shuffled
% rank data
% - pval: the non-parametric p-value

% To see the code in action for random data, run:
%{
dum = repmat([1:nreg],[nruns,1]);
ranks = zeros(size(dum));
for j = 1:size(dum,1)
    ranks(j,:) = dum(j,randperm(nreg));
end
[dZ,pval] = rankClustering(ranks, 500);
%}

function [dZ,pval] = rankClustering(ranks, nperm)

if nargin < 2
    error('Not enough input arguments')
end

nruns = size(ranks,1);
nreg = size(ranks,2);

%% compute rank histograms

ls = 1:nreg; % ranks
hdata = zeros(nreg); % region * rank
for s = 1:nreg
    hdata(s,:) = hist(ranks(:,s),ls);
end

%% euclidean distance between all data points, for each region

d = rankDistance(hdata,nruns);

%% expected euclidean distance - generated from shuffled data

df_store = zeros(nperm,nreg);

for np = 1:nperm
    % create fake data set
    dum = repmat([1:nreg],[nruns,1]);
    % shuffle
    ranksf = zeros(size(dum));
    for j = 1:size(dum,1)
        ranksf(j,:) = dum(j,randperm(nreg));
    end
    
    hdataf = zeros(size(hdata));
    for i = 1:size(dum,2)
        hdataf(i,:) = hist(ranksf(:,i),ls);
    end
    
    % euclidean distance of fake data
    df_store(np,:) = rankDistance(hdataf,nruns);
end

% z-score
dZ = (d-nanmean(df_store,1)') ./ nanstd(df_store,0,1)';

% pval
[~,id] = sort([d'; df_store],1);
loc = zeros(nreg,1);
for r = 1:nreg
    loc(r) = find(id(:,r)==1);
end
pval = loc/(nperm+1);

end


function dist = rankDistance(hdata,nruns)

dist = zeros(size(hdata,1),1);
for r1 = 1:size(hdata,2) % rank1
    for r2 = r1+1:size(hdata,2) % rank2
        dist = dist + abs(r1-r2)*hdata(:,r1).*hdata(:,r2);
    end
end
dist = dist/nruns^2;
end