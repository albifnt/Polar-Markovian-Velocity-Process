function A = matrix_A( Grid, Fluid, Perm, bCond )
nx = Grid.nx;
ny = Grid.ny;
del_x = Grid.del_x;
del_y = Grid.del_y;
del_z = Grid.del_z;
mu = Fluid.mu;
k_avg = Perm.k_avg;
A=zeros(nx*ny,nx*ny);
for q=1:nx*ny
    [i, j]=get_ij(q,nx,ny);
    if(i>1)
        k_=2/(1/k_avg(j,i-1)+1/k_avg(j,i)); % harmonic more accurate
        A(q,q-1)=k_*del_y*del_z/del_x/mu;
        A(q,q)=A(q,q)-A(q,q-1);
    end
    if(i<nx)
        k_=2/(1/k_avg(j,i)+1/k_avg(j,i+1));
        A(q,q+1)=k_*del_y*del_z/del_x/mu;
        A(q,q)=A(q,q)-A(q,q+1);
    end
    if(j>1)
        k_=2/(1/k_avg(j,i)+1/k_avg(j-1,i));
        A(q,q-nx)=k_*del_x*del_z/del_y/mu;
        A(q,q)=A(q,q)-A(q,q-nx);
    end
    if(j<ny)
        k_=2/(1/k_avg(j+1,i)+1/k_avg(j,i));
        A(q,q+nx)=k_*del_x*del_z/del_y/mu;
        A(q,q)=A(q,q)-A(q,q+nx);
    end
end
A = matrix_A_boundary(A, Grid, Fluid, Perm, bCond);
end