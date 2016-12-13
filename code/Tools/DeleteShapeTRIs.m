function shape_reduced = DeleteShapeTRIs(shape,tris_to_keep)
%%
shape_reduced = shape;

shape_reduced.TRIV = shape.TRIV(tris_to_keep==1 ,:);

vertices_to_keep = unique(shape_reduced.TRIV(:));
shape_reduced.X = shape.X(vertices_to_keep);
shape_reduced.Y = shape.Y(vertices_to_keep);
shape_reduced.Z = shape.Z(vertices_to_keep);

[isMem, shape_reduced.TRIV] = ismember(shape_reduced.TRIV,vertices_to_keep);
if ~all(isMem(:))
    error('bad triangulation')
end