% Given a flag matrix of a shape builds its laplacian operator with
% dirichlet boundary conditions
% L - Constrained laplacian. Rows / Cols - the location for each member

function [L rows cols dims new_index] = neumann3D (flag, d)

[row col dim] = size(flag);
number_of_pixels = sum(sum(sum(flag)));
rows = zeros(number_of_pixels,1);
cols = zeros(number_of_pixels,1);
L = sparse(number_of_pixels,number_of_pixels);
count = 1;
new_index = create_index(flag);


for i=1:row
    for j=1:col
        for k=1:dim
            
            if flag(i,j,k) < 1
                continue;
            end
            rows(count) = i;
            cols(count) = j;
            dims(count) = k;
            
            L(count,new_index(i,j,k)) = 6;
            
            if i==1 || flag(i-1,j,k)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i-1,j,k)) = -1;
            end
            if i==row || flag(i+1,j,k)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i+1,j,k)) = -1;
            end
            
            if j==1 || flag(i,j-1,k)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i,j-1,k)) = -1;
            end
            if j==col || flag(i,j+1,k)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i,j+1,k)) = -1;
            end
            
            
            if k==1 || flag(i,j,k-1)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i,j,k-1)) = -1;
            end
            if k==dim || flag(i,j,k+1)==0
                L(count,new_index(i,j,k)) = L(count,new_index(i,j,k)) - 1;
            else
                L(count,new_index(i,j,k+1)) = -1;
            end
            
            
            count = count + 1;
            
        end
    end
end

L = L / d^2;

end



% create new index for each shape member
function index = create_index(flag)
[row col dim] = size(flag);
index = zeros(row,col,dim);
count = 1;
for i=1:row
    for j=1:col
        for k=1:dim
            if flag(i,j,k) > 0
                index(i,j,k) = count;
                count = count + 1;
            end
        end
    end
end
end


