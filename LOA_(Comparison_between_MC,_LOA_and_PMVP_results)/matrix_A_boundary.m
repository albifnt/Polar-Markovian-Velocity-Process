function A = matrix_A_boundary(A, Grid, Fluid, Perm, bCond )
nx = Grid.nx;
ny = Grid.ny;
del_x = Grid.del_x;
del_y = Grid.del_y;
del_z = Grid.del_z;
mu = Fluid.mu;
k_avg = Perm.k_avg;
if(bCond.cond(1)==1)
    i=1;
    for j=1:ny
        k_ = k_avg(j,i);
        q = i+(j-1)*nx;
        A(q,q)=A(q,q)-2*k_*del_y*del_z/del_x/mu;
    end
end
if(bCond.cond(2)==1)
    i=nx;
    for j=1:ny
        k_ = k_avg(j,i);
        q = i+(j-1)*nx;
        A(q,q)=A(q,q)-2*k_*del_y*del_z/del_x/mu;
    end
end
if(bCond.cond(3)==1)
    j=1;
    for i=1:nx
        k_ = k_avg(j,i);
        q = i+(j-1)*nx;
        A(q,q)=A(q,q)-2*k_*del_x*del_z/del_y/mu;
    end
end
if(bCond.cond(4)==1)
    j=ny;
    for i=1:nx
        k_ = k_avg(j,i);
        q = i+(j-1)*nx;
        A(q,q)=A(q,q)-2*k_*del_x*del_z/del_y/mu;
    end
end
end