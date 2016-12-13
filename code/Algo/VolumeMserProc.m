function [mserRes shape] = VolumeMserProc(shape_file_name,params)
%% init
tStart = now;

shape_folder = params.data_path;
shape_full_path = fullfile(shape_folder,shape_file_name); 
assert(exist(shape_full_path,'file')>0,'shape file "%s" not found',shape_file_name)
[~,shape_name] = fileparts(shape_full_path);

%% try to load MSERs from mat-file
p             = fullfile(shape_folder ,'mser_res',params.benchName);
mser_res_file = fullfile(p,[shape_name '.mser.mat']);
if params.load_msers && exist(mser_res_file,'file')
    VerboseDisp('loading MSER result for shape %s',shape_name)
    mserRes = load(mser_res_file);
    if nargout > 1
        shape = load(shape_full_path);
    end
    return
elseif ~isdir(p),mkdir(p)
end
clear p

%% load shape
VerboseDisp('loading shape %s',shape_name)
shape = load(shape_full_path);
nVertices = numel(shape.index);

vertexVol = shape.d^3 * ones(nVertices,1);
shape.volume = vertexVol;

%% get laplacian
laplacian_filename = sprintf('%s.laplacian.mat',shape_name);
laplacian_fullpath = fullfile(shape_folder ,'laplacian',laplacian_filename);
lb = GetShapeLB(laplacian_fullpath,shape,params.max_num_evals,params.load_laplacian);

ADJ = lb.L < 0;
adjLut = Adj2Lut(ADJ);

%% ComponentTree
VerboseDisp('building ComponentTree')

if params.is_ew
    ewFun = GetEdgeWeightFun(params.ew_distFun,shape,lb.evals,lb.evecs,params.timeVals,lb.siHksDesc,lb.siHksNorm);
    [ M, TT, area, val, idx ] = AgglomerativeTree(ADJ , vertexVol(:), ewFun ,params.ew_clusterType); %#ok<ASGLU>
else
    vw = GetVertexWeights(params.vw_fun,lb.evals,lb.evecs,params.timeVals,lb.siHksDesc,lb.siHksNorm);
    if params.vw_flip_sign,vw = -vw;end % find MSER "minus"
    [ M, TT, area, val, idx ] = ComponentTree( adjLut, vertexVol(:), vw(:) );%#ok<ASGLU>
end

shape.val = val;

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
    mserVertexArea = vertexVol( msers(iMser).isInMser );
    
    msers(iMser).instability    = mserProps.instability(   iMser);
    msers(iMser).val            = mserProps.val(           iMser);
    %msers(iMser).area           = mserProps.area(          iMser);
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
            lb.evals,...
            lb.evecs(isInMser(:,iMser),:),...
            mserVertexArea,...
            params.desc_vs);
    end
    
end


%%
VerboseDisp('Finalizing...')
save(mser_res_file,...
    'shape',...
    'msers')

mserRes.msers = msers;
mserRes.shape = shape;

VerboseDisp('finished shape %s, found %d MSERs, took %s \n',....
    shape_name,nMsers,datestr(now-tStart,13) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lb = GetShapeLB(laplacian_fullpath,shape,max_num_evals,isLoad)
%%
if isLoad && exist(laplacian_fullpath,'file')
    VerboseDisp('loading laplacian...')
    lb = load(laplacian_fullpath);
    
else
    VerboseDisp('calculating laplacian...')
    [lb.L, ~, ~, ~, IDX] = neumann3D (shape.G, shape.d);
    perm = IDX(shape.G>0);
    lb.L = lb.L(perm,perm);
    
    [lb.evecs lb.evals] = prepare_evec(lb.L, max_num_evals);
    
    p = fileparts(laplacian_fullpath);
	
    MakeDir(p)

    [lb.siHksDesc ,lb.siHksNorm] = Eigen2SIHKS(lb.evals,lb.evecs);
    
    save(laplacian_fullpath,'-struct','lb','L','evecs','evals','siHksDesc','siHksNorm')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shape_reduced = DeleteShapeVertices(shape,vertices_to_keep)
%% find TRIs who's 3 vertices are"to keep"
isMem = ismember(shape.TRIV,vertices_to_keep);
tris_to_keep = all(isMem,2);
shape_reduced = DeleteShapeTRIs(shape,tris_to_keep);


