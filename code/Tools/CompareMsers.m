function [mats res1 res2] = CompareMsers(msers1,msers2,vertexArea1,lut2to1)
% msers1        - region struct of 1st shape (master)
% msers2        - region struct of 2nd shape (slave)
% vertexArea1   - area of vertices in master shape
% lut2to1       - look-up table for the location of slave's vertices in master

isOut3 = nargout > 1;

mats.instability1 = [msers1.instability];
mats.instability2 = [msers2.instability];
if isOut3
    res1.instability = mats.instability1;
    res2.instability = mats.instability2;
end

nMsers1 = numel(msers1);
nMsers2 = numel(msers2);
intersectionRatio = zeros(nMsers1,nMsers2);

% pre-calc mser1's area...
mserArea1 = zeros(nMsers1,1);
for iMser1 = 1:nMsers1
    mserArea1(iMser1) = sum(vertexArea1(msers1(iMser1).isInMser));
end


for iMser2 = 1:nMsers2
    mser2Vertices = lut2to1(msers2(iMser2).isInMser);% = mser2 vertices in shape1 indexes
    mserArea2 = sum(vertexArea1(mser2Vertices));
    
    for iMser1 = 1:nMsers1
        %
        isInMser1 = msers1(iMser1).isInMser(mser2Vertices); % flag if mser2 vertex is also in mser1
        if ~any(isInMser1); continue; end
        intersectingVertices = mser2Vertices(isInMser1);
        intersectionArea     = sum(vertexArea1(intersectingVertices));
        
        % intersectionArea cannot be bigger than either areas (happens due to "lut2to1" non-injective-ness)
        intersectionArea     = min(intersectionArea,mserArea1(iMser1));
        intersectionArea     = min(intersectionArea,mserArea2);
        
        intersectionRatio(iMser1,iMser2) = intersectionArea / ...
            (mserArea1(iMser1) + mserArea2 - intersectionArea);
        
    end
end
mats.intersectionRatio = intersectionRatio;


if ~nMsers1 || ~nMsers2
    mats.desc_dist = [];
    if  isOut3
        res1.maxIntersection   = zeros(size(res1.instability));
        res2.maxIntersection   = zeros(size(res2.instability));
        res1.matchIntersection = zeros(size(res1.instability));
        res2.matchIntersection = zeros(size(res2.instability));
    end
    return 
    
elseif isOut3
    % keep each MSER's max intersection
    res1.maxIntersection = max(intersectionRatio,[],2)';
    res2.maxIntersection = max(intersectionRatio,[],1);
    
end

%% keep each MSER's NN intersection
msers1_desc = [msers1.desc];
msers2_desc = [msers2.desc];
descTypes = fieldnames(msers1_desc);

for iType = 1:numel(descTypes)
    currType = descTypes{iType};
    
    curr_desc1 = [msers1_desc.(currType)];
    curr_desc2 = [msers2_desc.(currType)];
    
    curr_desc_size1 = sum(curr_desc1.^2);
    curr_desc_size2 = sum(curr_desc2.^2);
    
    % find L2 distance by cosine-law (faster than ANN over all)
    desc_dist = bsxfun(@plus,curr_desc_size1',curr_desc_size2) - ...
        2 * curr_desc1' * curr_desc2;
    
    mats.desc_dist.(currType) = desc_dist;
    
    if isOut3
        [dummy,nnIdx] = min(desc_dist,[],2); %#ok<ASGLU>
        idxs = sub2ind([nMsers1,nMsers2],1:nMsers1,nnIdx');
        res1.matchIntersection.(currType) = intersectionRatio(idxs);
        
        [dummy,nnIdx] = min(desc_dist,[],1); %#ok<ASGLU>
        idxs = sub2ind([nMsers1,nMsers2],nnIdx,1:nMsers2);
        res2.matchIntersection.(currType) = intersectionRatio(idxs);
    end
    
end
%%

