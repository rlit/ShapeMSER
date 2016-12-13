function shape = loadoff(filename)

shape = [];

f = fopen(filename, 'rt');

n = '';
while isempty(n)
    fgetl(f);
    n = sscanf(fgetl(f), '%d %d %d');
end
    
nv = n(1);
nt = n(2);
data = fscanf(f, '%f');

shape.TRIV = reshape(data(end-3*nt+1:end), [3 nt])';


data = data(1:end-3*nt);
data = reshape(data, [length(data)/nv nv]);

shape.X = data(1,:)';
shape.Y = data(2,:)';
shape.Z = data(3,:)';

fclose(f);
