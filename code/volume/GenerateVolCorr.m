function volCorr = GenerateVolCorr(xform,null,corr)
%%
xform.idxsSurface = GetVolSurface(xform.G);
null.idxsSurface  = GetVolSurface( null.G);

xform.idxsInterior = setdiff(xform.index ,xform.idxsSurface);
null.idxsInterior  = setdiff( null.index , null.idxsSurface);

%% Interior - Surface Corr 
[~,xform.i2s] = InteriorToSurfaceCorr(xform.G);
persistent G s2i
if ~isequal(G,null.G)
    G = null.G;
    s2i = InteriorToSurfaceCorr(G);
end
null.s2i = s2i;

%% create surface correspondance
surfaceCorr = GetSurfaceCorr(null,xform,corr);

%% handle parts of null surface that are not on the lut "null.s2i" 
surfaceCorr = FixSurfaceCorr(surfaceCorr,null);

%% create interior correspondance
interiorCorr = GetInteriorCorr(surfaceCorr,xform,null);

%% unite surface & interior 
volCorr = zeros(size(xform.index));

[~, locs] = ismember(xform.idxsSurface,xform.index);
volCorr(locs) = surfaceCorr;

[~, locs] = ismember(xform.idxsInterior,xform.index);
volCorr(locs) = interiorCorr;

volCorr = struct(...
    'xformIdx',num2cell(xform.index),...
    'nullIdx',num2cell(volCorr));

function surfaceCorr = GetSurfaceCorr(null,xform,corr)
%%
[xform.idxsSurface,surface_xform] = GetVolSurface(xform.G,xform.g);

tree = ann('init',[xform.X xform.Y xform.Z ]');
closestVertex = ann('search', tree, surface_xform' , 1)';
ann('deinit', tree);

if isnumeric(corr) && isvector(corr)
    closestNullIdx = corr(closestVertex);
    closestNullLoc = [null.X(closestNullIdx) null.Y(closestNullIdx) null.Z(closestNullIdx) ];
else
    % lut = lookup table from xform vertex to null xyz location
    if isnumeric(corr) && size(corr,2)==6
        lut = CreateLutFromCorr(null,xform,corr); 
    else
        lut = corr;
    end
    assert(size(lut,2)==3)
    closestNullLoc = lut(closestVertex,:);
end

[null.idxsSurface,surface_null] = GetVolSurface(null.G,null.g);

tree = ann('init',surface_null');
surfaceCorr = ann('search', tree, closestNullLoc' , 1)';
ann('deinit', tree);

surfaceCorr = null.idxsSurface(surfaceCorr);

function surfaceCorr = FixSurfaceCorr(surfaceCorr,null)
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


for iPoint = 1:numel(notOnS2i)
    toBeChanged = surfaceCorr == notOnS2i(iPoint);
    if ~any(toBeChanged),continue,end
    
    shortestPath  = shortestpath(T,tmpSrc(iPoint,:));
    if isempty(shortestPath)
        newI = ann('search', tree, tmpSrc(iPoint,:)' , 1)';
        warning('ugly patch used')
    else
        newI = ann('search', tree, shortestPath(end,:)' , 1)';
    end
    
    surfaceCorr(toBeChanged) = null.s2i(newI).index;
end
ann('deinit', tree);

function interiorCorr = GetInteriorCorr(surfaceCorr,xform,null)
%%
nInteriorVoxels = numel(xform.idxsInterior);

xfrom2Surface = zeros(size(xform.G));
xfrom2Surface(xform.idxsSurface) = 1:numel(xform.idxsSurface);
nulltoS2i = zeros(size(null.G));
nulltoS2i([null.s2i.index]) = 1:numel(null.s2i);

avgDistDiff = 0;
maxDistDiff = 0;
interiorCorr = zeros(nInteriorVoxels,1);
for ii = 1:nInteriorVoxels
    xformDist  = xform.i2s(ii).distance;
    xform_Idx  = xform.i2s(ii).lut;
    xform_sIdx = xfrom2Surface(xform_Idx);
    
    null_sIdx = surfaceCorr(xform_sIdx);
    null_sIdx = nulltoS2i(null_sIdx);
    
    nullDist = null.s2i(null_sIdx).distance;
    distDiff = abs(xformDist - nullDist);
    [~,minLoc] = min(distDiff);
    
    avgDistDiff =     avgDistDiff + distDiff(minLoc);
    maxDistDiff = max(maxDistDiff , distDiff(minLoc));
    
    interiorCorr(ii) = null.s2i(null_sIdx).lut(minLoc);
end
avgDistDiff = avgDistDiff / numel(xform.idxsInterior);

