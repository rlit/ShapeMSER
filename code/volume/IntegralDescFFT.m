function desc = IntegralDescFFT(volIndicator,radii)

% TODO - check smoothing filter size
% dilateRadius = 3;
% dilateFilter = ones(dilateRadius,dilateRadius,dilateRadius) / dilateRadius^3;
% shapeG_dilated = double(convn(volIndicator,dilateFilter,'same'));
shapeG_dilated = volIndicator;

G = fftn(shapeG_dilated);
[s1,s2,s3] = size(G);
[idx2,idx1,idx3] = meshgrid(1:s2,1:s1,1:s3) ;



% desc = zeros(nnz(volIndicator),numel(radii));
desc = cell(numel(radii),1);
for i = 1:numel(radii)
    %t=tic;
    curRadius = radii(i);
    
    % create sphere filter
    curFilter = (...
        (idx1 - s1/2).^2 + ...
        (idx2 - s2/2).^2 + ...
        (idx3 - s3/2).^2) < curRadius^2;
    curFilter = curFilter / nnz(curFilter);
    F = abs(fftn(curFilter)); % the "abs" is to prevent any shift of the filter

    % apply filter on volmue matrix
    filterRes = ifftn(G.*F);
    %desc{i} = filterRes(volIndicator~=0);
    desc{i} = filterRes;
    
    
    %disp([curRadius toc(t)])          
end

% %%
% isosurface(curFilter,.5/5000)
% axis equal
% axis([1 s2 1 s1 1 s3])
% %%
% sum(abs(real(F(:))))
