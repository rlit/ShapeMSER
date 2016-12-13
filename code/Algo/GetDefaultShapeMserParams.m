function params = GetDefaultShapeMserParams(params)

if ~exist('params','var') || ~isstruct(params),params = struct;end

%% main function flags
if ~isfield(params,'showPlots'),   params.showPlots    = 0;     end
if ~isfield(params,'preCalcMsers'),params.preCalcMsers = 1;     end

if ~isfield(params,'benchName'),params.benchName = 'default';   end


%% pre-component tree params - mainly LBO's
% pre-calc usage flags
if ~isfield(params,'load_laplacian'),params.load_laplacian = 1; end
if ~isfield(params,'load_msers'),    params.load_msers     = 0; end

if ~isfield(params,'maxVertices'),  params.maxVertices   = 1e4; end
if ~isfield(params,'min_comp_size'),params.min_comp_size = .05; end
if ~isfield(params,'max_boundry_dist'),params.max_boundry_dist = 5; end

if ~isfield(params,'removeZeroEval'),params.removeZeroEval = 1; end
if ~isfield(params,'max_num_evals'),params.max_num_evals = 200; end
if ~isfield(params,'cut_regions'),  params.cut_regions = 0;     end



%% component tree params
if ~isfield(params,'is_ew'), params.is_ew  = 1; end % <-- vw / ew flag
if ~isfield(params,'timeVals'),params.timeVals = [1024 2048];end

% ------------ vw - vertex-weight ----------------------------------------
if ~isfield(params,'vw_fun'),
    %for all options see file "GetVertexWeights.m"
    params.vw_fun  = 'heat_kernel';
end
if ~isfield(params,'vw_flip_sign'),params.vw_flip_sign = 0;     end

% ------------ ew - edge-weight -------------------------------------------
if ~isfield(params,'ew_clusterType'),
    %for all options see file "AgglomerativeTree.m"
    params.ew_clusterType  = 'min';
end
if ~isfield(params,'ew_distFun'), 
    %for all options see file "GetEdgeWeightFun.m"
    params.ew_distFun  = 'inverse_HKS';
end


%% mser detection params
if ~isfield(params,'mser_filters')
    params.mser_filters = struct(...        
        'MinPeakClimb', 1,... <-- disabled if ge(1)
        'KeepStable', 0,...
        'MinArea', 0.01, ...
        'MaxArea', 0.75, ...
        'MaxScore', inf,...
        'MinDiversity', 0.2 );
end

%% descriptors 
if ~isfield(params,'calc_desc'),params.calc_desc = true; end %vocabulary size
if ~isfield(params,'desc_vs'),params.desc_vs = 6:2:12; end %vocabulary size