function [G_ index X Y Z]= sample_boundary(G)

G_=0*G;

[row col dim] = size(G);
index = zeros(row,col,dim);
count = 1; index=[];
X=[];Y=[];Z=[];
 for i=1:row
  for j=1:col
      for k=1:dim
          
          
          if G(i,j,k)==0
              continue;
          end
         
               
          if i==1 || G(i-1,j,k)==0  || i==row || G(i+1,j,k)==0 ||  j==1 || G(i,j-1,k)==0 ...
                      ||   j==col || G(i,j+1,k)==0 ||  k==1 || G(i,j,k-1)==0  ...
                         || k==dim || G(i,j,k+1)==0      
         
                      
                     if G(i,j,k) > 0
                        index=[index; count];
                        G_(i,j,k) = 1;
                        X=[X;i];Y=[Y;j]; Z=[Z;k];
                       
                     end
          end
          
           count = count + 1;
          
      end
  end
 end
end

