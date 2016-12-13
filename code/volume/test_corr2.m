
null  = load('0001.null.0.vol.mat');
shapeName =  '0001.affine.5.vol.mat';
xform = load(shapeName);
corrLut  = load(strrep(shapeName,'mat','corr.mat'));
corrLut    = [corrLut.volCorr.nullIdx]';
[~,corrLut] = ismember(corrLut,find(null.G));

lb = load(strrep(shapeName,'mat','laplacian.mat'));
%%

t = 4000;
K = (lb.evecs.^2)*exp(-2*lb.evals(:)*t);

% Build adjacency
ADJ = double(lb.L<0);
[ ~, TT, area, val, idx] = ComponentTree( ADJ, ones(size(ADJ,1),1)/size(ADJ,1), K(:) );
[mserIdxs] = FindMSERs(TT, val, area, struct('MinArea', 0.01, 'MaxArea', 0.75, 'MinDiversity', 0.2) );
nMsers = numel(mserIdxs);

%%
close all;

figure;
subplot(122)
h = plot3(xform.points(:,1),xform.points(:,2),xform.points(:,3), '.k');
set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
axis image;
hold on

subplot(121)
h = plot3(null.points(:,1),null.points(:,2),null.points(:,3), '.k');
set(h,'Marker','o', 'MarkerSize', 1, 'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
axis image;
hold on

for i = 1:nMsers,
    isInMser = FloodFill(ADJ, K(:), idx(mserIdxs(i)));
    
    subplot(122)
    h_(1) = plot3(xform.points(isInMser,1),xform.points(isInMser,2),xform.points(isInMser,3), '.k');
    set(h_(1),'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
    title(sprintf('region #%d ',i))
    
    
    isInNull = false(size(null.index));
    isInNull(unique( corrLut(isInMser))) = true;
    nnz(isInNull)
    
    subplot(121)
    h_(2) = plot3(null.points(isInNull,1),null.points(isInNull,2),null.points(isInNull,3), '.r');
    set(h_(2),'Marker','o', 'MarkerSize', 3, 'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0]);
        
    pause;
    delete(h_)
end


