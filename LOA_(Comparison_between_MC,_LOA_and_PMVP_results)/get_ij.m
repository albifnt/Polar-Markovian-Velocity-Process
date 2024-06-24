function [i, j]=get_ij(q,nx,ny)
i=rem(q,nx);
j=fix(q/nx)+1;
if(i==0)
    i=nx;
    j=j-1;
end
end