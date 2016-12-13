function [s2i i2s] = InteriorToSurfaceCorr(G)
%% get shape surface
idxsSurface = GetVolSurface(G);
[subsSurface(:,1) subsSurface(:,2) subsSurface(:,3) ] = ind2sub(size(G),idxsSurface);

%% run FMM from surfave inwards
F = G;
F(F==0) = 1e-8;%dont march through voxels out of the shape
T = msfm(F, subsSurface');

%%
idxsInterior = setdiff(find(G),idxsSurface);
[subsInterior(:,1) subsInterior(:,2) subsInterior(:,3) ] = ind2sub(size(T),idxsInterior);

%%
closestPoints = zeros(size(subsInterior));
if matlabpool('size') ~= 0
    FindClosestSurface_PAR()
else
    FindClosestSurface()
end

%% convert balck to indexes
tree = ann('init',subsSurface');
closestSurface = ann('search', tree, closestPoints' , 1);
ann('deinit', tree);


%%
i2s = struct(...
    'index',   num2cell(idxsInterior),...
    'lut',     num2cell(idxsSurface(closestSurface)),...
    'distance',num2cell(T(idxsInterior)));

%%
s2i = struct('index',   num2cell(idxsSurface));
toDelete = false(size(s2i));
for ii = 1:numel(s2i)
    curI = closestSurface == ii;
    s2i(ii).lut = idxsInterior(curI);
    if isempty(s2i(ii).lut)
        toDelete(ii) = true;
        continue
    end
    s2i(ii).distance = T(s2i(ii).lut);
end
s2i = s2i(~toDelete);

%%
% figure
% h = plot3(subsSurface(:,1),subsSurface(:,2),subsSurface(:,3), '.k');
% set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% hold on
% h = plot3(subsSurface(~toDelete,1),subsSurface(~toDelete,2),subsSurface(~toDelete,3), '.r');
% set(h,'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
% axis image;
% hold off;


%%
    function FindClosestSurface()
        
        iTree = ann('init',subsInterior');
        
        isDone = false(size(idxsInterior))';        
        [dummy,perm] = sort(T(idxsInterior),'descend');
        for iPoint = 1:size(subsInterior,1)
            curPoint = perm(iPoint);
            
            if isDone(curPoint), continue,end            
            shortestPath  = shortestpath(T,subsInterior(curPoint,:));
            closestPoints(curPoint,:) = shortestPath(end,:);
            isDone(curPoint) = true;
            
            [lineIdx,lineD] = ann('search', iTree, shortestPath' , 1);
            lineIdx = lineIdx(lineD<0.5 & ~isDone(lineIdx));
            lineIdx = unique(lineIdx);

            closestPoints(lineIdx,:) = repmat(shortestPath(end,:),numel(lineIdx),1);
            isDone(lineIdx) = true;
            
        end
        ann('deinit', iTree);        

    end

%%
    function FindClosestSurface_PAR()
        parfor iPoint = 1:size(subsInterior,1)
            shortestPath  = shortestpath(T,subsInterior(iPoint,:));
            closestPoints(iPoint,:) = shortestPath(end,:);
        end
    end
end