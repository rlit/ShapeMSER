if ~exist('FindMSERs','file')
    p = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(p,'InitPath.m'))
end

% Load shape
if ~exist('shape','var'),
    clear shape;
    shape = loadoff('0001.null.0.off');
%     shape = loadoff('0002.isometry.3.off');
    
    d = 2;
    [G d g points index] = mesh2volume_D(shape,d);
    [L rows cols dims IDX] = neumann3D (G, d);
    perm = IDX(G>0);
    L = L(perm,perm);
    
    [evecs evals] = prepare_evec(L, 150);
end

%%
% close all;
t = 1000;
K = (evecs.^2)*exp(-2*evals(:)*t);
% 
% K_ = zeros(size(G));
% K_(G>0) = K; %(index);
% 
% figure
% VAL = smooth3(K_,'gaussian',5);
% [X,Y,Z] = meshgrid(g{2},g{1},g{3});
% h = slice(X,Y,Z,log10(VAL+1),[0 0],[42 -30],[73 73]);
% set(gca,'clim',prctile(VAL(:),[1 99]))
% set(h,'linestyle','none')
% axis image;
% % K = VAL(G>0);

%% Build adjacency
ADJ = double(L<0);
[ M, TT, area, val, idx, root ] = ComponentTree( ADJ, ones(size(ADJ,1),1)/size(ADJ,1), K(:) );
[mserIdxs, mserProps] = FindMSERs(TT, val, area, struct('MinArea', 0.01, 'MaxArea', 0.75, 'MinDiversity', 0.2) );
nMsers = numel(mserIdxs);

%%
figure;
h = plot3(points(:,1),points(:,2),points(:,3), '.k');
set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
axis image;
%MASK = zeros(size(G));
for i = 1:nMsers,
    isInMser = FloodFill(ADJ, K(:), idx(mserIdxs(i)));

    i1 = find(isInMser == 1);
    hold on
    h = plot3(points(i1,1),points(i1,2),points(i1,3), '.r');
    set(h,'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
    hold off;    
    title(sprintf('region #%d ',i))
    %saveas(gcf,int2str(i+30),'png')
    pause(1);
    delete(h)
end



