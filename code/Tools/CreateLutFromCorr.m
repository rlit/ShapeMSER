function lut = CreateLutFromCorr(null_shape,trans_shape,corr_data)
% lut = lookup table from xform vertex to null euclidean location
%%
% [dummy perm] = sort(corr_data(:,5));
% corr_data = corr_data(perm,:);

%%
nCorr = size(corr_data,1);
lut = zeros(nCorr,3);
for iCorr = 1:nCorr
    
    curr_null_tri         = corr_data(iCorr,1);
    curr_null_barycentric = corr_data(iCorr,2:4)';
    curr_null_vertices = null_shape.TRIV(curr_null_tri,:);
    
    curr_null_euclidean = sum([
        null_shape.X(curr_null_vertices) .* curr_null_barycentric ...
        null_shape.Y(curr_null_vertices) .* curr_null_barycentric ...
        null_shape.Z(curr_null_vertices) .* curr_null_barycentric]);
    
    %%
    
    curr_trans_tri = corr_data(iCorr,5);
    curr_trans_tri_vertex = corr_data(iCorr,6:8) == 1;
    if ~any(curr_trans_tri_vertex), error('barycentric coordinates are not on a vertex'),end
    
    curr_trans_vertex = trans_shape.TRIV(curr_trans_tri,curr_trans_tri_vertex);
    if any(lut(curr_trans_vertex,:)), error('attempt to update LUT in a non-empty location'),end
    
    lut(curr_trans_vertex,:) = curr_null_euclidean;
end

