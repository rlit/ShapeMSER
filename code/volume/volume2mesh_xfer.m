function s = volume2mesh_xfer(points, shape, v)

tree = ann('init', points');
[i,d] = ann('search', tree, [shape.X(:) shape.Y(:) shape.Z(:)]', 1, 'eps', 0);
ann('deinit', tree);
ann('close');

s = v(i,1);
%trisurf(shape.TRIV, shape.X, shape.Y, shape.Z, s); axis image; axis off; shading interp;