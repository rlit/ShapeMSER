function desc = IntegralDesc(volIndicator,scales)

% TODO - check smoothing filter size
dilateFilter = ones(3,3,3) / 5^3;
shapeG_dilated = double(convn(volIndicator,dilateFilter,'same'));

% desc = zeros(nnz(volIndicator),numel(scales));
desc = cell(numel(scales),1);
for i = 1:numel(scales)
    t=tic;
    curRadius = scales(i);
    
    % create sphere filter
    curVec = -curRadius:curRadius;
    [sphereX,sphereY,sphereZ] = meshgrid(curVec,curVec,curVec) ;
    sphereFilter = (...
        (sphereX).^2 + ...
        (sphereY).^2 + ...
        (sphereZ).^2) < curRadius^2;
    sphereFilter = sphereFilter / nnz(sphereFilter);
    
    
    % apply filter on volmue matrix
    filterRes = convn(shapeG_dilated,sphereFilter,'same');
    %desc(:,i) = filterRes(volIndicator~=0);
    desc{i} = filterRes;
    
    disp([curRadius toc(t)])
          
end

