function shrec2010_bench(params)
%% set params
if ~exist('params','var') || ~isstruct(params)
    params = GetDefaultBenchParams();
else
    params = GetDefaultBenchParams(params);
end

nShapes = 3;
nStrengh = 5;

transformationClasses = GetShrec2010ClassNames();
nClasses = numel(transformationClasses);

%% get data - assuming all needed files are in ..\data
if isempty(params.data_path)
    params.data_path = GetShrecFolderLocation('SHREC10');
end

plot_path = fullfile(params.data_path,'figures',params.benchName);
MakeDir(plot_path)

%% pre-calc MSERs on all shapes first (Parallelize-able)
if params.preCalcMsers && matlabpool('size') > 0
    PreCalcMsers(params)
    params.load_msers = true;
end

%%
showPlots = params.showPlots;
benchRes = cell2struct(cell(nClasses,nShapes),transformationClasses,1);
for iShape = 1:nShapes
    % read null shape
    nullShapeName = GetShrecShapeFileName(iShape,'null',0);
    [nullMserRes null_shape] = FileMserProc([nullShapeName '.off'],params);
    
    % points used to build kd-tree: from null shape euclidean to null vertex #
    shape_points = [...
        nullMserRes.shape.X ...
        nullMserRes.shape.Y ...
        nullMserRes.shape.Z ]';
    
    % load look-up table for bilateral symmetry
    symData = load([nullShapeName '.sym.mat']);
    nullSymLut = symData.isym(nullMserRes.shape.origIdx);
    null_tree = ann('init',shape_points);
    nullSymLut = ann('search', null_tree, [...
        null_shape.X(nullSymLut) ...
        null_shape.Y(nullSymLut) ...
        null_shape.Z(nullSymLut) ]' , 1)';
    ann('deinit', null_tree);
    
    % read transformed shapes
    shapeRes = cell(nClasses,1);
    for iClass = 1:nClasses
        currClass = transformationClasses{iClass};

        if showPlots
            f = figure(iClass + iShape * nClasses);
            clf(f)
            colormap([0.6 0.6 0.6]);
            subplot(2, 3, 1);cla reset;
            PlotShapeMsers(nullMserRes)
            title(['shape #' int2str(iShape)]);
        end
        
        for iStrengh = 1:nStrengh
            currShapeName = GetShrecShapeFileName(iShape,currClass,iStrengh);
            [currMserRes curr_shape] = FileMserProc([currShapeName '.off'],params);
            
            %% create shape-to-null LUT:
            shape_full_path = fullfile(params.data_path,[currShapeName '.off']);
            % 1) from curr shape vertices' index to null shape vertices' euclidean
            lut = CreateShapeLut(shape_full_path,curr_shape,null_shape);
            % 2) reduce lut to fit reduced transformed shape
            lut = lut(currMserRes.shape.origIdx,:);
            % 3) convert lut: from transformed shape index to null's
            null_tree = ann('init',shape_points);
            lut = ann('search', null_tree, lut', 1)';
            ann('deinit', null_tree);
            
            %% comapre every possible MSER pair between shapes
            classRes(iStrengh,1) = GetMserComparison(...
                nullMserRes.msers, ...
                currMserRes.msers, ...
                nullMserRes.shape.area, ...
                lut,nullSymLut(lut)); %#ok<AGROW>
            
            %% plot
            if showPlots
                subplot(2, 3, iStrengh + 1,'parent',f);
                PlotShapeMsers(currMserRes)
                title([currClass ' #' int2str(iStrengh)]);
            end
        end % end strengh loop
        shapeRes{iClass} = classRes;
        %% finallize class display
        if showPlots
            % to set colors right
            allAxes = findobj('type','axes','-and','parent',f);
            set(allAxes,'CLim',[0 size(colormap,1)-1])
            
            figName = [int2str(iShape) '_' currClass];
            figPath = fullfile(plot_path,figName);
            saveas(f,[figPath,'.jpg'])
            saveas(f,[figPath,'.fig'])
            
            if ishandle(f),close(f),end
        end
    end % end transform class loop
    % finallize shape
    benchRes(iShape) = cell2struct(shapeRes,transformationClasses,1);
    
end % end shape loop

%%
ann('close');

res_path = fullfile(params.data_path,'benchRes');

MakeDir(res_path)
save(fullfile(res_path,['benchRes_' params.benchName '.mat']),'benchRes','params')

% plot result
Plot_Bench_Results(benchRes,params.benchName,plot_path)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PreCalcMsers(params)
params.load_msers = false;
shapeFiles = dir(fullfile(params.data_path,'*.off'));
shapeFiles = {shapeFiles.name};
parfor iShape = 1:numel(shapeFiles)
    shapeName = shapeFiles{iShape};
    [dummy dummy] = FileMserProc(shapeName,params); %#ok<NASGU>
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotShapeMsers(mserRes,iMser)
if exist('iMser','var')
    mserRes.msers = mserRes.msers(iMser);
end
%% prep data
nMesrsOnVertex = sum(horzcat(mserRes.msers.isInMser),2);
maxVertexMsers = max(nMesrsOnVertex);
%% plot
if maxVertexMsers > 0
    shape_ = TessellateShape( mserRes.shape, nMesrsOnVertex );
else
    shape_ = mserRes.shape;
    shape_.tri_labels = nMesrsOnVertex;
end
PlotShape(shape_,shape_.tri_labels-1)

view([0 20]);
light;
lighting phong;
camlight head;

if size(colormap,1) < maxVertexMsers + 1
    colormap([0.6 0.6 0.6; jet(maxVertexMsers)]);
end
caxis([0 size(colormap,1)-1]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compareRes = GetMserComparison(nullMsers,currMsers, nullVertexArea, lut,lutSym)
% this function alters the outer variable "classRes"
% the reason why this funtion is nested & not appended in the end
% of the m-file is for cases that there are no msers detected
%% regular lut
[compareRes.mats compareRes.null compareRes.xform] = CompareMsers(...
    nullMsers, ...
    currMsers, ...
    nullVertexArea, ...
    lut);


%% symmetry - compare again, with the symmetic lut

[compareResSym.mats compareResSym.null compareResSym.xform  ] = CompareMsers(...
    nullMsers, ...
    currMsers, ...
    nullVertexArea, ...
    lutSym);

%% compare the two previous scores - take smaller distances & bigger ovelaps
compareRes.mats.intersectionRatio = max(compareResSym.mats.intersectionRatio, compareRes.mats.intersectionRatio);
compareRes.null.maxIntersection  = max(compareResSym.null.maxIntersection,  compareRes.null.maxIntersection);
compareRes.xform.maxIntersection = max(compareResSym.xform.maxIntersection, compareRes.xform.maxIntersection);

if ~isfield(compareRes.mats,'desc_dist') || isempty(compareRes.mats.desc_dist)
    return
end

descTypes = fieldnames(compareRes.mats.desc_dist);
for iType = 1:numel(descTypes)
    currType = descTypes{iType};
    compareRes.mats.desc_dist.(         currType) = min(compareResSym.mats.desc_dist.(currType),          compareRes.mats.desc_dist.(currType));
    compareRes.null.matchIntersection.( currType) = max(compareResSym.null.matchIntersection.(currType),  compareRes.null.matchIntersection.(currType));
    compareRes.xform.matchIntersection.(currType) = max(compareResSym.xform.matchIntersection.(currType), compareRes.xform.matchIntersection.(currType));    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
