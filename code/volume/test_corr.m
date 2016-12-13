if ~exist('FindMSERs','file')
    p = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(p,'InitPath.m'))
end

% Load shapes
if ~exist('null','var'),
    %%
    tic
    null  = loadoff('0001.null.0.off');
    xform = loadoff('0001.isometry.3.off');
    corr  = load('0001.isometry.3.corr');
    
    %npt   = 75;
    d = 2;
    [ null.G  null.d  null.g  null.points  null.index] = mesh2volume_D(null, d);
    [xform.G xform.d xform.g xform.points xform.index] = mesh2volume_D(xform,d);
    toc
    
end

%% create surface correspondance
[xform.idxsSurface,surface_xform] = GetVolSurface(xform.G,xform.g);
% % h = plot3(xform.points(:,1),xform.points(:,2),xform.points(:,3), '.k');
% % set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% % hold on
% % h = plot3(surface_xform(:,1),surface_xform(:,2),surface_xform(:,3), '.r');
% % set(h,'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
% % axis image;
% % hold off;

tree = ann('init',[xform.X xform.Y xform.Z ]');
closestVertex = ann('search', tree, surface_xform' , 1)';
ann('deinit', tree);

lut = CreateLutFromCorr(null,xform,corr); % lut = lookup table from xform vertex to null euclidean location
closestNullLoc = lut(closestVertex,:);

[null.idxsSurface,surface_null] = GetVolSurface(null.G,null.g);

tree = ann('init',surface_null');
surfaceCorrLut = ann('search', tree, closestNullLoc' , 1)';
ann('deinit', tree);


%%
% clf
% surface_null_offset = bsxfun(@plus,surface_null,[110 70 30]);
% h = plot3(surface_null_offset(:,1),surface_null_offset(:,2),surface_null_offset(:,3), '.k');
% set(h,'Marker','o', 'MarkerSize', 2, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% hold on
% h = plot3(surface_xform(:,1),surface_xform(:,2),surface_xform(:,3), '.r');
% set(h,'Marker','o', 'MarkerSize', 2, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
% axis image;
% 
% currI = 0;
% step = 100;
% while 1
%     currI = mod(currI,step)+1;
%     xformI = currI:step:numel(surfaceCorrLut);
%     nullI = surfaceCorrLut(xformI);
%     h = plot3(...
%         [surface_xform(xformI,1) surface_null_offset(nullI,1)]',...
%         [surface_xform(xformI,2) surface_null_offset(nullI,2)]',...
%         [surface_xform(xformI,3) surface_null_offset(nullI,3)]', 'g');
%     pause(5)
%     delete(h)
% end

%% Interior To Surface Corr
disp('InteriorToSurfaceCorr - null')
tic
[null.s2i null.i2s] = InteriorToSurfaceCorr(null.G);
toc
disp('InteriorToSurfaceCorr - xform')
tic
[xform.s2i xform.i2s] = InteriorToSurfaceCorr(xform.G);
toc

%% 

xform.idxsInterior = setdiff(xform.index ,xform.idxsSurface);
null.idxsInterior  = setdiff( null.index , null.idxsSurface);

%% handle parts of null surface that are not on the lut "null.s2i" 

tmpSrc = [];
[tmpSrc(:,1),tmpSrc(:,2),tmpSrc(:,3)] = ind2sub(size(null.G),[null.s2i.index]);
tree = ann('init',tmpSrc');

nullSurfaceIm = ones(size(null.G)) * 1e-8;
nullSurfaceIm(null.idxsSurface) = 1;

T = msfm(nullSurfaceIm, tmpSrc');


notOnS2i = setdiff(null.idxsSurface,[null.s2i.index]);
tmpSrc = [];
[tmpSrc(:,1),tmpSrc(:,2),tmpSrc(:,3)] = ind2sub(size(null.G),notOnS2i);

surfaceCorrLut_ = null.idxsSurface(surfaceCorrLut);

for iPoint = 1:numel(notOnS2i)
    toBeChanged = surfaceCorrLut_ == notOnS2i(iPoint);
    if ~any(toBeChanged),continue,end
    
    shortestPath  = shortestpath(T,tmpSrc(iPoint,:));
    newI = ann('search', tree, shortestPath(end,:)' , 1)';
    
    surfaceCorrLut_(toBeChanged) = null.s2i(newI).index;
end
ann('deinit', tree);

%%
finalLut = zeros(size(xform.index));
[~, locs] =ismember(xform.idxsSurface,xform.index);
finalLut(locs) = null.idxsSurface(surfaceCorrLut);

nulltoS2i = zeros(size(null.G));
nulltoS2i([null.s2i.index]) = 1:numel(null.s2i);

%%
[~, locs] =ismember(xform.idxsInterior,xform.index);
xfrom2Surface = zeros(size(xform.G));
xfrom2Surface(xform.idxsSurface) = 1:numel(xform.idxsSurface);

avgDistDiff = 0;
maxDistDiff = 0;
for ii = 1:numel(xform.idxsInterior)
    xformDist  = xform.i2s(ii).distance;
    xform_Idx  = xform.i2s(ii).lut;
    xform_sIdx = xfrom2Surface(xform_Idx);
    
    null_sIdx = surfaceCorrLut_(xform_sIdx);
    null_sIdx = nulltoS2i(null_sIdx);
    
    nullDist = null.s2i(null_sIdx).distance;
    distDiff = abs(xformDist - nullDist);
    [~,minLoc] = min(distDiff);
    
    avgDistDiff =     avgDistDiff + distDiff(minLoc);
    maxDistDiff = max(maxDistDiff , distDiff(minLoc));
    
    finalLut(locs(ii)) = null.s2i(null_sIdx).lut(minLoc);
end
avgDistDiff =     avgDistDiff / numel(xform.idxsInterior);
%%


% surfaceDist = VolSurfaceDist(null.G);
% for v = 1:max(surfaceDist(:))
%     
%     idxs = find(surfaceDist == v);
%     [I1,I2,I3] = ind2sub(size(null.G),idxs);
%     
%     I1 = null.g{1}(I1);
%     I2 = null.g{2}(I2);
%     I3 = null.g{3}(I3);
%     points = [I1(:) I2(:) I3(:)];
%     
%     surface_offset = bsxfun(@plus,points,[0 70*v 0]);
%     h = plot3(surface_offset(:,1),surface_offset(:,2),surface_offset(:,3), '.r');
%     set(h,'Marker','o', 'MarkerSize', 2, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
%     axis image;
%     hold on
% end

