function v = baricentric2vertex(shape,t,u)

[m,i] = max(u,[],2);
v = zeros(length(t),1);

for k = 1:length(t)
    v(k) = shape.TRIV(t(k),i(k));
end

