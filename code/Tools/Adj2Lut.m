function [lut edges] = Adj2Lut(adjacencyMtx)
nVertices = length(adjacencyMtx);

[r c] = find(tril(adjacencyMtx)); %<- tril = assuming non-directed graph
nEdges = length(r);

lut = cell(nVertices,1);

% % pre-allocate?
% nNeighbors   = accumarray([r;c],ones(1,nEdges*2));
% maxNeighbors = max(nNeighbors);

for iEdge = 1:nEdges
    currR = r(iEdge);
    currC = c(iEdge);
    lut{currR} = [lut{currR} currC];
    lut{currC} = [lut{currC} currR];
end

if nargout > 1
    edges = [r c];
end
