function [desc,norm] = Eigen2SIHKS(evals,evecs,timeSamples)

if ~exist('timeSamples','var'),
    timeSamples = 2.^(1: 1/16 : 25-1/16);
end

% avgDesc is the average descriptor on a representative dataset
persistent avgDesc %<-- to prevent loading more than once
if isempty(avgDesc)
    loadRes = load('SIHKS_avg_desc.mat');
    avgDesc = loadRes.MN;
end
%%
desc = scalespace(timeSamples, evals, evecs,1/2);
lastNonZero = find(all(desc),1, 'last');
desc = desc(:,1:lastNonZero);

% remove sacling by derivative of a log
desc = diff(log(desc),1,2);

% remove shift by taking abs of ftt
desc = abs(fft(desc,[],2));

% subsract avgDesc
desc = bsxfun(@minus,desc,avgDesc(1:lastNonZero-1));

if nargout > 1
    norm = sum(desc.^2,2);
end

desc = desc(:,1:6);



