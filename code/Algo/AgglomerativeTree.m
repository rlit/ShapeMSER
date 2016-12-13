function [ parentIdxs, TT, nodeArea, f, idx, root ] = AgglomerativeTree( adjacency, a, distFun ,clusterType)

if exist('clusterType','var')
    switch clusterType
        case 'min'
            % no need for joinThFun
        case 'avg'
            joinThFun = @mean;
        case 'med'
            joinThFun = @median;
        case '2nd'
            joinThFun = @(x)x(min(numel(x),2));
        case '20%'
            joinThFun = @(x)prctile(x,20);            
        otherwise
            clusterType = 'min';
    end
else
    clusterType = 'min';
end

%%
a = a / sum(a);
nVertices = length(a);
isSetInterior = false(nVertices,1);

% find all edges & their distances
[node1 node2] = find(tril(adjacency)); %<- tril = assuming non-directed graph
nEdges = length(node1);

% calc edge_weights in steps to prevent "Out of memory"
edge_weights = zeros(nEdges,1);
stepSize = 1e4;
for iStep = 0:floor(nEdges/stepSize)
    currRange =    (1 + iStep  * stepSize) : ...
        min(nEdges,(1 + iStep) * stepSize);
    edge_weights(currRange) = distFun(node1(currRange),node2(currRange));
end, clear stepSize iStep currRange

[edge_weights,perm] = sort(edge_weights);
node1 = node1(perm);
node2 = node2(perm);

edgeList = [node1 node2];
edgeLut  = Edge2Lut(edgeList,nVertices);

disjointSets = unionfind('init', 2*nVertices);
setRoot   = 1:2*nVertices;

parentIdxs   = zeros(2*nVertices,1);
childrenIdxs = zeros(2*nVertices,2);
nChildren    = ones(2*nVertices,1);

nodeArea = [a(:) ; zeros(nVertices,1)];
f = zeros(2*nVertices,1);

isEdgeJoined = false(size(edge_weights));
for iEdge = 1:nEdges
    if isEdgeJoined(iEdge),continue,end
    
    % find sets # of nodes
    curSet1 = unionfind('find', disjointSets, node1(iEdge));
    curSet2 = unionfind('find', disjointSets, node2(iEdge));
    curDist = edge_weights(iEdge);
    %% check if nodes are alredy connected
    if curSet1 == curSet2
        continue
    end
    %% find set root
    curRoot1 = setRoot(curSet1);
    curRoot2 = setRoot(curSet2);
    
    % sanity check - set-roots can't have parents (yet)
    if parentIdxs(curRoot1) || parentIdxs(curRoot2)
        error('sanity check failed')
    end
    
    %% comapre iEdge to all other edges connecting the 2 sets
    if ~strcmp(clusterType,'min')
        % find all edges connecting the 2 sets
        if nChildren(curRoot1) > nChildren(curRoot2)
            connEdges = FindConnectingEdges(curRoot2,curSet1);
        else
            connEdges = FindConnectingEdges(curRoot1,curSet2);
        end
        connDists = sort(edge_weights(connEdges));
        
        % (if not pass threshold - skip iEdge)
        distTh = joinThFun(connDists);
        if isempty(distTh),distTh=0;end
        
        if curDist < distTh
            continue
        end
    end
    
    
    %% perform "join" on the 2 sets
    %add a parent to the 2 nodes
    nVertices = nVertices + 1;
    parentIdxs(curRoot1) = nVertices;
    parentIdxs(curRoot2) = nVertices;
    childrenIdxs(nVertices,:) = [curRoot1 curRoot2];
    
    % sum areas
    nodeArea(nVertices)  = nodeArea(curRoot1)  + nodeArea(curRoot2);
    nChildren(nVertices) = nChildren(curRoot1) + nChildren(curRoot2);
    f(nVertices)         = curDist;
    
    %mark the 2 roots & the parent as "joined"
    joinedSet = MergeSets(disjointSets,  curSet1, curSet2);
    joinedSet = MergeSets(disjointSets,  joinedSet, nVertices);
    
    % remember the set's root
    setRoot(joinedSet) = nVertices;
    isEdgeJoined(iEdge) = true;
    
end
%% free memory
parentIdxs(nVertices+1:end) = [];
nodeArea(  nVertices+1:end) = [];
f(         nVertices+1:end) = [];
unionfind('deinit', disjointSets);
clear disjointSets

%% finalize
idx = 1:nVertices;
root = (parentIdxs == 0); % NOTE - # of roots equals the # of connected component
TT = sparse(...
    parentIdxs(~root),...
    idx(~root),...
    ones(nVertices-nnz(root),1), ...
    nVertices,nVertices);
root = find(root);

%% NESTED FUNCTIONS
    function connEdges = FindConnectingEdges(root,targetSet)
        connEdges = [];
        rootSet   = unionfind('find', disjointSets, root);
        leaves    = GetNodeLeaves(root);
        
        for iLeaf = 1:numel(leaves)
            currLeaf = leaves(iLeaf);
            if isSetInterior(currLeaf),continue,end
            leafEdges = edgeLut{currLeaf};
            
            isAllEdgesJoined = true;
            for iLeafEdge = 1:numel(leafEdges)
                currLeafEdge = leafEdges(iLeafEdge);
                if isEdgeJoined(currLeafEdge),continue,end
                
                otherNode = edgeList(currLeafEdge,:);
                otherNode = otherNode(otherNode ~= currLeaf);
                otherSet  = unionfind('find', disjointSets, otherNode);
                
                if otherSet == rootSet
                    isEdgeJoined(currLeafEdge) = true; continue
                end
                isAllEdgesJoined = false;
                
                if otherSet == targetSet
                    connEdges = [connEdges currLeafEdge]; %#ok<AGROW>
                end
                
            end % leaf-edges loop
            isSetInterior(currLeaf) = isAllEdgesJoined;
        end %leaves loop
        connEdges = unique(connEdges);
    end

    function leaves = GetNodeLeaves(nodes)
        leaves  = [];
        while ~isempty(nodes)
            currChildren = childrenIdxs(nodes(:),:);
            isLeaf = all(currChildren==0,2);
            leaves = [leaves ; nodes(isLeaf)]; %#ok<AGROW>
            nodes  = currChildren(currChildren > 0);
        end
        
    end

end

function joinedSet = MergeSets(qnode,  n1, n2)
unionfind('union', qnode, n1, n2);
joinedSet = unionfind('find', qnode, n1);
end

function edgeLut  = Edge2Lut(edgeList,nVertices)
edgeLut = cell(nVertices,1);
for iEdge = 1:size(edgeList,1)
    edgeNode1 = edgeList(iEdge,1);
    edgeNode2 = edgeList(iEdge,2);
    edgeLut{edgeNode1} = [edgeLut{edgeNode1} iEdge];
    edgeLut{edgeNode2} = [edgeLut{edgeNode2} iEdge];
end
end
