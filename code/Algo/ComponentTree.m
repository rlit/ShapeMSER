function [ M, TT, area, f, idx, root ] = ComponentTree( ADJ, a, f )
if iscell(ADJ)
    adjLut = ADJ;
else
    adjLut = Adj2Lut(ADJ);
end
a = a / sum(a);
N = length(a);

qnode = unionfind('init', N);
qtree = unionfind('init', N);
lowestNode = (1:N);
level      = f(:);
highest    = f(:);
area       = a(:);

T = sparse([],[],[],N,N,N);

[f_,idx] = sort(f, 'descend');
isProcessed = false(N,1);
for k = 1:length(idx),
    p = idx(k);
    curTree = unionfind('find', qtree, p);
    curNode = unionfind('find', qnode, lowestNode(curTree));
    
    % idx1 = processed neighbors
    currNeighbors = adjLut{p};
    idx1 = currNeighbors(isProcessed(currNeighbors));

%    fprintf(1, 'Step %d\n', k);
    
    for j = 1:length(idx1),
        q = idx1(j);
        adjTree = unionfind('find', qtree, q);
        adjNode = unionfind('find', qnode, lowestNode(adjTree));
        
        if curNode == adjNode, continue; end
        
        if level(curNode) == level(adjNode)
            % Merge
%            fprintf(1, '   Merge:    %d (%g)   %d (%g)\n', curNode, f(curNode), adjNode, f(adjNode) );
            [curNode, area, highest] = MergeNodes(qnode, area, highest, adjNode, curNode);
        else
            % Add child 
%            fprintf(1, '   Child:    %d (%g) -> %d (%g)\n', curNode, f(curNode), adjNode, f(adjNode) );
            T(curNode,adjNode) = 1;
            area(curNode) = area(curNode) + area(adjNode);
            highest(curNode) = max( highest(curNode), highest(adjNode) );
        end
        
%        fprintf(1, '   Mtree:    %d %d\n', curTree, adjTree );
        unionfind('union', qtree, adjTree, curTree);
        
        curTree = unionfind('find', qtree, curTree);
        T(lowestNode(curTree),lowestNode(adjTree)) = 1;
        lowestNode(curTree) = curNode;
                
    end
    
    isProcessed(p) = true;
end

root = lowestNode(unionfind('find', qtree, unionfind('find', qnode, 1)));
M = zeros(N,1);
for p = 1:N,
    M(p) = unionfind('find', qnode, p);
end

[Mu,idx,idx1] = unique(M, 'first');
TT = sparse(length(Mu),length(Mu));

[i,j] = find(T);
ii = unique([idx1(i),idx1(j)],'rows');
i = ii(:,1); j = ii(:,2);
k = find(i ~= j);
i = i(k); j = j(k);
TT = sparse(i,j,ones(length(i),1), length(Mu),length(Mu));

unionfind('deinit', qnode);
unionfind('deinit', qtree);
clear qnode qtree;

root = idx1(root);
f = f(idx);
area = area(idx);

