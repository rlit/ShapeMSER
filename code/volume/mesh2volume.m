function [G d grid points idx] = mesh2volume(shape,N,show)

show = 0;
if nargin == 3
    show = 1;
end

mX = min(shape.X);  MX = max(shape.X);
mY = min(shape.Y);  MY = max(shape.Y);
mZ = min(shape.Z);  MZ = max(shape.Z);

d = max([ (MX-mX)/N (MY-mY)/N (MZ-mZ)/N]);
n = [(MX-mX)/d (MY-mY)/d (MZ-mZ)/d];
grid{1} = linspace(mX,MX,n(1));
grid{2} = linspace(mY,MY,n(2));
grid{3} = linspace(mZ,MZ,n(3));

G = voxelize( grid{1} , grid{2}, grid{3} ,  shape.X(shape.TRIV)',shape.Y(shape.TRIV)',shape.Z(shape.TRIV)');
G = SmoothVolumeSurface(G);

%grid{2} = grid{2}(30:end-31);
%grid{3} = grid{3}(23:end-23);
%G = G(:,30:end-31,23:end-23);

[I,J,K] = ndgrid(grid{1}, grid{2} , grid{3});

idx = find(G==1);
points = [I(idx) J(idx) K(idx)];

if show == 1
    plot3(I(idx),J(idx),K(idx),'k.'), axis image
end

