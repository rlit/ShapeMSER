% Merge nodes n1, n2
function [tmpNode, area, highest] = MergeNodes(qnode, area, highest, n1, n2)

unionfind('union', qnode, n1, n2);
tmpNode = unionfind('find', qnode, n1);
if tmpNode == n2,
    tmpNode2 = n1;
else
    tmpNode2 = n2;
end
area(tmpNode)    = area(tmpNode) + area(tmpNode2);
highest(tmpNode) = max( highest(tmpNode), highest(tmpNode2) );
