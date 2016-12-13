function branchRes = Climb(parentsIdxs, node)

% branchRes = mxClimb(parentsIdxs, node)';

%% old matlab code that did the same as mxClimb
%% stupid little code to prevent calling "zeros" mony times
persistent branchAllocated allocSize
if ~isequal(allocSize,numel(parentsIdxs))
    allocSize = numel(parentsIdxs);
    branchAllocated = zeros(allocSize,1);
end
%% main code of "climb"
for idx = 1:numel(parentsIdxs)
    branchAllocated(idx) = node;
    parent = parentsIdxs(node);
    if ~parent ,break; end
    node = parent;
end

% reverse order & remove zeros
branchRes = flipud(branchAllocated(1:idx));
