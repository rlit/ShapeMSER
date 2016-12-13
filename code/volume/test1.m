if 1,
addpath ../iso2mesh/
addpath ../common/

shape = loadoff('0001.null.0_.off');

bb0 = min([shape.X shape.Y shape.Z],[],1);
bb1 = max([shape.X shape.Y shape.Z],[],1);

volume = [];
[vert,volume.elem,volume.TETV] = ...
    surf2mesh([shape.X shape.Y shape.Z], shape.TRIV, bb0,bb1,0.05,25);
volume.X2 = vert(:,1);
volume.Y2 = vert(:,2);
volume.Z2 = vert(:,3);

[A B] = build_fem_3D(volume);

num_evecs = 200;
[evecs, evals] = eigs(A, B, num_evecs, 'sm'); %0, struct('disp', 0));
evals = diag(evals);

[evals, idx] = sort(evals);
evecs = evecs(:,idx);
end

% Translate on surface
tree = ann('init', [volume.X2(:) volume.Y2(:) volume.Z2(:)]');
[idx, dist] = ann('search', tree, [shape.X(:) shape.Y(:) shape.Z(:)]', 1, 'eps', 0);
ann('deinit', tree);
ann('close');


[z,i] = max(volume.Z2);
Kxx = [];
Kxy = [];
for t=[50 100 200 500], 
    Kxy(:,end+1) = bsxfun(@times, evecs(i,:), evecs) * exp(-evals*t);
    Kxx(:,end+1) = (evecs.^2) * exp(-evals*t);
end

for k=1:size(Kxx,2),
    %figure(1);
    %subplot(2,2,k);
    %trisurf(volume.TETV(:,1:3), volume.X2, volume.Y2, volume.Z2, Kxx(:,k)); axis image;

    figure(1);
    subplot(2,2,k);
    trisurf(shape.TRIV,shape.X,shape.Y,shape.Z, Kxx(idx,k)); axis image; shading interp;

    figure(2);
    subplot(2,2,k);
    trisurf(shape.TRIV,shape.X,shape.Y,shape.Z, Kxy(idx,k)); axis image; shading interp;

end
