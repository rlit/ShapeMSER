function [surfaceDist,idxs] = VolSurfaceDist(volMat)
%%
f(:,:,1) = [0 0 0 ; 0 1   0 ; 0 0 0];
f(:,:,2) = [0 1 0 ; 1 0.5 1 ; 0 1 0];
f(:,:,3) = [0 0 0 ; 0 1   0 ; 0 0 0];
%%
convRes = convn(volMat,f,'same');
isSurface = rem(convRes,1)==0.5 & convRes<6;
%%
surfaceDist = double(isSurface);
f(2,2,2) = 0;
layerNumber = 2;

while 1
    convRes = convn(surfaceDist,f,'same');
    innerLayer = volMat>0 & convRes>0 & surfaceDist==0;
    if ~any(innerLayer(:))
        break
    end
    
    surfaceDist(innerLayer) = layerNumber;
    layerNumber = layerNumber + 1;

end


