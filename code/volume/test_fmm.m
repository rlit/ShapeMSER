shape = loadoff('0001.isometry.5.off');
%%
d = 1;
[G d g points index] = mesh2volume_D(shape,d);
%%
tic
[s2i i2s] = InteriorToSurfaceCorr(G);
toc
%%
[sIdxs,sPoints] = GetVolSurface(G,g);
[sSubs(:,1) sSubs(:,2) sSubs(:,3) ] = ind2sub(size(G),sIdxs);

%%
F = G;
F(F==0) = 1e-8;
tic
T = msfm(F, sSubs');
toc
% T(T>1e7) = NaN;

%%
% [maxLoc] = find(T==max(T(T<1e7)));
% T(maxLoc)
% [startPoints(1,:) startPoints(2,:) startPoints(3,:) ] = ind2sub(size(T),maxLoc);

%%
[iIdxs orig] = setdiff(find(G),sIdxs);
% iPoints = points(orig,:);
[startPoints(:,1) startPoints(:,2) startPoints(:,3) ] = ind2sub(size(T),iIdxs);

%%
tic
closestPoints = zeros(size(startPoints));
parfor iPoint = 1:size(startPoints,1)
    ShortestLine  = shortestpath(T,startPoints(iPoint,:));
    closestPoints(iPoint,:) = ShortestLine(end,:);
end
toc
%%
tree = ann('init',sSubs');
[closestVertex d] = ann('search', tree, closestPoints' , 1);
ann('deinit', tree);
%%
hist(double(closestVertex),2000)
% nnz(closestVertex==426)

%%
figure(1);clf
h = plot3(sSubs(:,1),sSubs(:,2),sSubs(:,3), '.k');
set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
axis image;


for iPoint = 1:50:size(startPoints,1)
    %               shortestpath(DistanceMap,StartPoint,SourcePoint,Stepsize,Method)
    ShortestLine  = shortestpath(T,startPoints(iPoint,:));
    
    tree = ann('init',sSubs');
    [closestVertex d] = ann('search', tree, ShortestLine(end,:)' , 1);
    ann('deinit', tree);

    hold on
    H(1) = plot3(sSubs(closestVertex,1),sSubs(closestVertex,2),sSubs(closestVertex,3), '.r');
    set(H(1),'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
    
    %H(2) = plot3(startPoints(iPoint,1), startPoints(iPoint,2), startPoints(iPoint,3), '.g');
    H(2) = plot3(ShortestLine(:,1),ShortestLine(:,2),ShortestLine(:,3), '.g');
    set(H(2),'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 1 0]);
    hold off;
    pause
    delete(H)
end

%%
p = [];
D = -0;
[p(:,1) p(:,2) p(:,3) ] = ind2sub(size(T),find(T>D & T<=D+1));

figure(1);clf
h = plot3(p(:,1) ,p(:,2), p(:,3), '.k');
set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
axis image;

