function params = GetDefaultBenchParams(params)
%% get algo prams
if ~exist('params','var') || ~isstruct(params)
    params = GetDefaultShapeMserParams();
else
    params = GetDefaultShapeMserParams(params);
end
%% add main function flags
if ~isfield(params,'showPlots'),   params.showPlots    = 0;     end
if ~isfield(params,'preCalcMsers'),params.preCalcMsers = 1;     end

if ~isfield(params,'benchName'),params.benchName = 'default';   end
if ~isfield(params,'data_path'),params.data_path = '';   end

