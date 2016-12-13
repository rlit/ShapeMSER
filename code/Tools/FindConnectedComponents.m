function [componentIdx,nComponents] = FindConnectedComponents(ADJ,min_comp_size)

if nargin < 2
    min_comp_size = 0;
end

nVertices = length(ADJ);

componentIdx = zeros(nVertices,1);
nComponents  = 0;
while ~all(componentIdx)
    % flood fill the 1st vertex with (componentIdx == 0)
    isInComponent = FloodFill(ADJ, ones(nVertices,1), find(~componentIdx,1));
    if (nnz(isInComponent)/nVertices > min_comp_size)
        nComponents = nComponents + 1;
        componentIdx(isInComponent) = nComponents;
    else
        componentIdx(isInComponent) = -1;
    end
end