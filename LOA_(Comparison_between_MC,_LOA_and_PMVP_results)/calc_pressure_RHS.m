function b = calc_pressure_RHS( Grid, Fluid, Perm, bCond )
nx = Grid.nx;
ny = Grid.ny;
del_x = Grid.del_x;
del_y = Grid.del_y;
del_z = Grid.del_z;
mu = Fluid.mu;
k_avg = Perm.k_avg;
b=sparse(nx*ny,1);
i=1;
for j=1:ny
    k_ = k_avg(j,i);
    q = i+(j-1)*nx;
    if(bCond.cond(1)==1)
        b(q) = b(q) - 2*bCond.value(1)*k_*del_y*del_z/del_x/mu;
    else
        b(q) = b(q) - bCond.value(1)*del_y*del_z;
    end
end
i=nx;
for j=1:ny
    k_ = k_avg(j,i);
    q = i+(j-1)*nx;
    if(bCond.cond(2)==1)
        b(q) = b(q) - 2*bCond.value(2)*k_*del_y*del_z/del_x/mu;
    else
        b(q) = b(q) - bCond.value(2)*del_y*del_z;
    end
end
j=1;
for i=1:nx
    k_ = k_avg(j,i);
    q = i+(j-1)*nx;
    if(bCond.cond(3)==1)
        b(q) = b(q) - 2*bCond.value(3)*k_*del_x*del_z/del_y/mu;
    else
        b(q) = b(q) - bCond.value(3)*del_x*del_z;
    end
end
j=ny;
for i=1:nx
    k_ = k_avg(j,i);
    q = i+(j-1)*nx;
    if(bCond.cond(4)==1)
        b(q) = b(q) - 2*bCond.value(4)*k_*del_x*del_z/del_y/mu;
    else
        b(q) = b(q) - bCond.value(4)*del_x*del_z;
    end
end
end