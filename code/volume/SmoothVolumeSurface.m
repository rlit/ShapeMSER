function volMat = SmoothVolumeSurface(volMat,minNeighbors)
if nargin < 2
    minNeighbors = 3;
end

%%
f(:,:,1) = [0 0 0 ; 0 1   0 ; 0 0 0];
f(:,:,2) = [0 1 0 ; 1 0.5 1 ; 0 1 0];
f(:,:,3) = [0 0 0 ; 0 1   0 ; 0 0 0];

%%
convRes = convn(volMat,f,'same');

toRemove = rem(convRes,1)==0.5 & convRes<minNeighbors;

volMat(toRemove) = 0;

