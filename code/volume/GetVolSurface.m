function [idxs,locations] = GetVolSurface(volMat,meshGrid)
%%
f(:,:,1) = [0 0 0 ; 0 1   0 ; 0 0 0];
f(:,:,2) = [0 1 0 ; 1 0.5 1 ; 0 1 0];
f(:,:,3) = [0 0 0 ; 0 1   0 ; 0 0 0];

%%
convRes = convn(volMat,f,'same');
isSurface = rem(convRes,1)==0.5 & convRes<6;

%%
idxs = find(isSurface);

if nargout > 1
    [I1,I2,I3] = ind2sub(size(volMat),idxs);
    locations = [I1(:) I2(:) I3(:)];
    
    if nargin > 1        
        I1 = meshGrid{1}(I1);
        I2 = meshGrid{2}(I2);
        I3 = meshGrid{3}(I3);
        locations = [I1(:) I2(:) I3(:)];
    end
    
end
