function b = calc_CYP_RHS( Grid, Fluid, Perm, bCond, i_Y, j_Y, p0)
nx = Grid.nx;
ny = Grid.ny;
nxperm = Grid.nxperm;
nyperm = Grid.nyperm;
del_x = Grid.del_x;
del_y = Grid.del_y;
del_z = Grid.del_z;
mu = Fluid.mu;
k_avg = Perm.k_avg;
CYY = Perm.CYY;
b=zeros(nx*ny,1);
for i=1:nx
    for j=1:ny
        q = i+(j-1)*nx;
        if(i > 1 && i < nx)
            J_n = (p0(q) - p0(q-1))/del_x;
            J_p = (p0(q+1) - p0(q))/del_x;
            CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i-1, 2*j, nxperm, nyperm);
            CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i+1, 2*j, nxperm, nyperm);
            k_p = 2/(1/k_avg(j,i+1)+1/k_avg(j,i));
            k_n = 2/(1/k_avg(j,i-1)+1/k_avg(j,i));
            b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_y*del_z/mu;
        end
        if(j > 1 && j < ny)
            J_n = (p0(q) - p0(q-nx))/del_y;
            J_p = (p0(q+nx) - p0(q))/del_y;
            CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j-1, nxperm, nyperm);
            CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j+1, nxperm, nyperm);
            k_p = 2/(1/k_avg(j+1,i)+1/k_avg(j,i));
            k_n = 2/(1/k_avg(j-1,i)+1/k_avg(j,i));
            b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_x*del_z/mu;
        end
    end
end
i=1;
for j=1:ny
    q = i+(j-1)*nx;
    if(bCond.cond(1) == 1)
        J_n = 2*(p0(q) - bCond.value(1))/del_x; 
    else
        J_n = - bCond.value(1)*mu/k_avg(j,i)/del_x;
    end
    J_p = (p0(q+1) - p0(q))/del_x;
    CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i-1, 2*j, nxperm, nyperm);
    CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i+1, 2*j, nxperm, nyperm);
    k_p = 2/(1/k_avg(j,i+1)+1/k_avg(j,i));
    k_n = k_avg(j,i);
    b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_y*del_z/mu;
end
i=nx;
for j=1:ny
    q = i+(j-1)*nx;
    J_n = (p0(q) - p0(q-1))/del_x;
    if(bCond.cond(2) == 1)
        J_p = 2*(bCond.value(2)-p0(q))/del_x;
    else
        J_p = bCond.value(2)*mu/k_avg(j,i)/del_x;
    end
    CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i-1, 2*j, nxperm, nyperm);
    CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i+1, 2*j, nxperm, nyperm);
    k_p = k_avg(j,i);
    k_n = 2/(1/k_avg(j,i-1)+1/k_avg(j,i));
    b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_y*del_z/mu;
end
j=1;
for i=1:nx
    q = i+(j-1)*nx;
    if(bCond.cond(3) == 1)
        J_n = 2*(p0(q) - bCond.value(3))/del_y;
    else
        J_n = - bCond.value(3)*mu/k_avg(j,i)/del_y;
    end
    J_p = (p0(q+nx) - p0(q))/del_y;
    CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j-1, nxperm, nyperm);
    CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j+1, nxperm, nyperm);
    k_p = 2/(1/k_avg(j+1,i)+1/k_avg(j,i));
    k_n = k_avg(j,i);
    b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_x*del_z/mu;
end
j=ny;
for i=1:nx
    q = i+(j-1)*nx;
    J_n = (p0(q) - p0(q-nx))/del_y;
    if(bCond.cond(4) == 1)
        J_p = 2*(bCond.value(4)-p0(q))/del_y;
    else
        J_p = bCond.value(4)*mu/k_avg(j,i)/del_y;
    end
    CYY_n = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j-1, nxperm, nyperm);
    CYY_p = get_CYY(CYY, i_Y, j_Y, 2*i, 2*j+1, nxperm, nyperm);
    k_p = k_avg(j,i);
    k_n = 2/(1/k_avg(j-1,i)+1/k_avg(j,i));
    b(q) = b(q) - (k_p*J_p*CYY_p - k_n*J_n*CYY_n)*del_x*del_z/mu;
end
end