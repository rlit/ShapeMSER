
p = GetDefaultBenchParams();
p.data_path = GetDataFolderLocation();
p.calc_desc = 0;

p.timeVals = 2048;

[mserRes shape_loaded] = FileMserProc('centaur4.mat',p);

%% plot result 
% prep data
nMesrsOnVertex = sum(horzcat(mserRes.msers.isInMser),2);
maxVertexMsers = max(nMesrsOnVertex);

% plot
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

colormap([0.6 0.6 0.6]);
if size(colormap,1) < maxVertexMsers + 1
    colormap([0.6 0.6 0.6; jet(maxVertexMsers)]);
end
caxis([0 size(colormap,1)-1]);