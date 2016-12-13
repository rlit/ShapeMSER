function [boundryVertices boundryEdges] = GetShapeBoundryVertices(shape)
allEdges = [...
    shape.TRIV(:,[1 2]) ; ...
    shape.TRIV(:,[1 3]) ; ...
    shape.TRIV(:,[2 3])];
allEdges = [min(allEdges,[],2) max(allEdges,[],2)];

allEdges = sortrows(allEdges);

isDiff = any(diff(allEdges),2);

boundryEdges = allEdges([1;isDiff] & [isDiff;1],:);
boundryVertices = unique(boundryEdges(:));

%boundryCoords = [shape.X(boundryVertices) shape.Y(boundryVertices) shape.Z(boundryVertices)];

