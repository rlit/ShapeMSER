function mserRes = ShapeMserProc(shape,params,shape_name)
%
%
%

nVertices = length(shape.X);

%% calc adjacency
ii = [shape.TRIV(:,[1 2]); shape.TRIV(:,[1 3]); shape.TRIV(:,[2 3])];
ii = unique([ii(:,[1 2]); ii(:,[2 1])], 'rows');
ADJ = sparse(ii(:,1), ii(:,2), ones(size(ii,1),1), nVertices, nVertices);
adjLut = Adj2Lut(ADJ);
clear ii


%% find connected components & run on each seperetally
[componentIdx,nComponents] = FindConnectedComponents(adjLut,params.min_comp_size);
shape.componentIdx = componentIdx;

if nComponents > 1
    %% this section is implemented but not tested. un comment to use at your own risk
    error ('shape conatains multiple components  - not supported yet')
%     VerboseDisp('Running on %d connected component(s)',nComponents)
%     
%     comp_msers = cell(nComponents,1);
%     shape.area = zeros(nVertices,1);
%     shape.val  = zeros(nVertices,1);
%     
%     for iComp = 1:nComponents
%         % create component shape
%         comp_name = sprintf('%s.%d',shape_name,iComp);
%         isInComp = (componentIdx == iComp);
%         comp_shape = DeleteShapeVertices(shape,find(isInComp));
%         
%         msers_res = ShapeMserProc(comp_shape,params,comp_name);
%         comp_msers{iComp} = msers_res.msers';
%         shape.val( isInComp) = mserRes.shape.val;
%         shape.area(isInComp) = mserRes.shape.area;
%         
%         VerboseDisp('Finised component %d/%d...',iComp,nComponents)
%     end
%     msers = [comp_msers{:}];

else
    %% get laplacian
    if params.load_laplacian && nargin > 2
        laplacian_filename = sprintf('%s.laplacian%d.mat',shape_name,nVertices);
        laplacian_fullpath = fullfile(params.data_path ,'laplacian',laplacian_filename);
    else
        laplacian_fullpath = '';
    end
    [evals,evecs,area,siHksDesc,siHksNorm] = GetShapeLB(...
        shape,...
        params.max_num_evals,...
        laplacian_fullpath);
    shape.area = area;
    
    %% ComponentTree
    VerboseDisp('building ComponentTree')
    
    if params.is_ew
        ewFun = GetEdgeWeightFun(params.ew_distFun,shape,evals,evecs,params.timeVals,siHksDesc,siHksNorm);
        [ M, TT, area, val, idx ] = AgglomerativeTree(ADJ , area(:), ewFun ,params.ew_clusterType); %#ok<ASGLU>
    else
        vw = GetVertexWeights(params.vw_fun,evals,evecs,params.timeVals,siHksDesc,siHksNorm);
        if params.vw_flip_sign,vw = -vw;end % find MSER "minus"
        [ M, TT, area, val, idx ] = ComponentTree( adjLut, area(:), vw(:) );%#ok<ASGLU>
    end
    
    shape.val = val(1:nVertices);
    
    %% Detect MSERs
    VerboseDisp('Detecting MSERs')
    [ mserIdxs, mserProps ] = FindMSERs(TT, val, area, params.mser_filters);
    nMsers = numel(mserIdxs);
    
    %% Flood fill to find MESR's vertices
    VerboseDisp('Flood Filling %d MSERs',nMsers)
    if params.is_ew
        isInMser = FloodFill(TT+TT',-val, idx(mserIdxs));
        isInMser = isInMser(1:nVertices,:);
    else
        isInMser = FloodFill(adjLut, vw(:), idx(mserIdxs));
    end
    
    
    %% Finalize
    VerboseDisp('Calc Descs...')
    msers = struct(...
        'isInMser',false(nVertices,1),...
        'area',false(nVertices,1),...
        'instability',0);
    msers = repmat(msers,nMsers,1);
    for iMser = 1:nMsers
        msers(iMser).isInMser = isInMser(:,iMser);
        mserVertexArea = shape.area( msers(iMser).isInMser );
        
        msers(iMser).instability    = mserProps.instability(   iMser);
        msers(iMser).val            = mserProps.val(           iMser);
        msers(iMser).peakClimbRatio = mserProps.peakClimbRatio(iMser);
        
        msers(iMser).area = sum(mserVertexArea);
        
        % % optional - "cut" mser & calc its laplacian
        %     mserShape = geodesic_contour(shape, ~msers(iMser).isInMser, 0.5);
        %     mserShape = DeleteShapeTRIs(mserShape,mserShape.tri_labels);
        %     nVertices = length(mserShape.X);
        %     [evecs, evals] = main_mshlp('cotangent', mserShape, min(nVertices,200));
        %     msers(iMser).evals = evals;
        %     msers(iMser).evecs = evecs;
        if params.calc_desc
            msers(iMser).desc = GetRegionDescs(...
                evals,...
                evecs(isInMser(:,iMser),:),...
                mserVertexArea,...
                params.desc_vs);
        end
        
    end
end

%%

mserRes.msers = msers;
mserRes.shape = shape;

% VerboseDisp('finished shape %s, found %d MSERs, took %s \n',....
%     shape_name,nMsers,datestr(now-tStart,13) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evals,evecs,area,siHksDesc,siHksNorm] = GetShapeLB(shape,max_num_evals,laplacian_fullpath)
%%

if exist(laplacian_fullpath,'file')
    VerboseDisp('loading laplacian...')
    load(laplacian_fullpath);
    
else
    VerboseDisp('calculating laplacian...')
    p = fileparts(laplacian_fullpath);
    MakeDir(p)
    nVertices = length(shape.X);
    [evecs, evals, W, A] = main_mshlp('cotangent', shape,  min(nVertices,max_num_evals)); %#ok<ASGLU>
    area = full(diag(A));
    %area = area/sum(area);
    
    [siHksDesc ,siHksNorm] = Eigen2SIHKS(evals,evecs);
    
    if isdir(fileparts(laplacian_fullpath))
        save(laplacian_fullpath,'evecs','evals','area','siHksDesc','siHksNorm')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function descStruct = GetRegionDescs(evals,evecs,vertex_weight,vocab_sizes)
regionArea = sum(vertex_weight);
%% SI-HKS
siHksDesc     = Eigen2SIHKS(evals,evecs);
for vs = vocab_sizes;
    vsStr = ['SI_vs_' int2str(vs)];
    trainRes = load(['vocab_' vsStr '.mat']);
    vqRes = SoftVQ(siHksDesc,...
        trainRes.vocab,...
        trainRes.sigma);
    
    mserSiDesc = vqRes * vertex_weight; % sum-weight by vertex vertex_weight
    descStruct.(vsStr) = mserSiDesc / regionArea;
end

mserSiDesc = siHksDesc' * vertex_weight;
descStruct.siDesc = mserSiDesc / regionArea;

%% regular HKS

for vs = vocab_sizes;
    vsStr = ['vs_' int2str(vs)];
    trainRes = load(['vocab_' vsStr '.mat']);
    hksDesc = Eigen2HKS(...
        evals,evecs,...
        trainRes.timeSamples);
    vqRes = SoftVQ(hksDesc,trainRes.vocab,trainRes.sigma);
    
    mserDesc = vqRes * vertex_weight; % sum-weight by vertex vertex_weight
    descStruct.(vsStr) = mserDesc / regionArea;
    
end

mserDesc = hksDesc' * vertex_weight;
descStruct.hksDesc = mserDesc / regionArea;
