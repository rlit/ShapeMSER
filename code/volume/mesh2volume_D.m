function [G d grid points idx] = mesh2volume_D(shape,d,show)


mX = min(shape.X);  MX = max(shape.X);
mY = min(shape.Y);  MY = max(shape.Y);
mZ = min(shape.Z);  MZ = max(shape.Z);

grid{1} = mX-d:d:MX+d;
grid{2} = mY-d:d:MY+d;
grid{3} = mZ-d:d:MZ+d;

G = voxelize( grid{1} , grid{2}, grid{3} ,  shape.X(shape.TRIV)',shape.Y(shape.TRIV)',shape.Z(shape.TRIV)');
G = SmoothVolumeSurface(G);

%grid{2} = grid{2}(30:end-31);
%grid{3} = grid{3}(23:end-23);
%G = G(:,30:end-31,23:end-23);

[I,J,K] = ndgrid(grid{1}, grid{2} , grid{3});

idx = find(G==1);
points = [I(idx) J(idx) K(idx)];

if nargin == 3
    plot3(I(idx),J(idx),K(idx),'k.'), axis image
end

