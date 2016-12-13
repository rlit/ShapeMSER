function X = barycentric2xyz(t,u,shape)

XX  = [shape.X shape.Y shape.Z];

tx = [XX(shape.TRIV(:,1),1) XX(shape.TRIV(:,2),1) XX(shape.TRIV(:,3),1)];
ty = [XX(shape.TRIV(:,1),2) XX(shape.TRIV(:,2),2) XX(shape.TRIV(:,3),2)];
tz = [XX(shape.TRIV(:,1),3) XX(shape.TRIV(:,2),3) XX(shape.TRIV(:,3),3)];

X = [dot(tx(t,:)',u')' dot(ty(t,:)',u')' dot(tz(t,:)',u')'];

