function shape_reduced = DeleteShapeVertices(shape,vertices_to_keep)
% find TRIs who's 3 vertices are "to keep"
isMem = ismember(shape.TRIV,vertices_to_keep);
tris_to_keep = all(isMem,2);
shape_reduced = DeleteShapeTRIs(shape,tris_to_keep);