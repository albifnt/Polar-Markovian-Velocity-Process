function vy_var = calc_vy_var( Grid, Fluid, Perm, bCond, p0, CYP, CPP )
nx = Grid.nx;
ny = Grid.ny;
nxperm = Grid.nxperm;
nyperm = Grid.nyperm;
del_y = Grid.del_y;
phi = Grid.poro;
mu = Fluid.mu;
k_avg = Perm.k_avg;
p0 = reshape(p0,nx,ny)';
vy_var = zeros((nx)*ny,1);
CYY = Perm.CYY;
for i=1:nx
    for j=1:ny
        q = i+(j-1)*(nx);
        if( j == 1 )
            k_ = k_avg(j, i);
            if bCond.cond(3) == 1 
                dpdy = 2*(p0(j,i) - bCond.value(3))/(del_y);
                vy_var(q) = vy_var(q) + (k_/mu/phi)^2 * dpdy^2 * get_CYY(CYY, 2*i, 2*j-1, 2*i, 2*j-1, nxperm, nyperm);
            else
                vy_var(q)=0;
            end
        elseif( j == ny ) 
            k_ = k_avg(j, i);
            if bCond.cond(4) == 1
                dpdy = 2*(bCond.value(4)-p0(j,i))/(del_y);
                vy_var(q) = vy_var(q) + (k_/mu/phi)^2*dpdy^2*get_CYY(CYY, 2*i, 2*j+1, 2*i, 2*j+1, nxperm, nyperm);
            else
                vy_var(q) = 0;
            end
        else
            k_ = k_avg(j, i);
            qz = i + (j-1) * nx; 
            qn = i + (j-2) * nx;
            qp = i + (j) * nx;
            dcpp_n = (CPP(qz,qn) - CPP(qn,qn))/del_y;
            dcpp_p = (CPP(qp,qp) - CPP(qz,qp))/del_y;
            ddcpp = (dcpp_p - dcpp_n)/del_y;
            vy_var(q) = vy_var(q) + (k_/mu/phi)^2 * ddcpp;
            dpdy = (p0(j+1,i) - p0(j-1,i))/(2*del_y);
            vy_var(q) = vy_var(q) + (k_/mu/phi)^2 * dpdy^2 * get_CYY(CYY, 2*i, 2*j, 2*i, 2*j, nxperm, nyperm);
            CYP_p = get_CYP(CYP, i, j+1, 2*i, 2*j, nx, ny, nxperm, nyperm);
            CYP_n = get_CYP(CYP, i, j-1, 2*i, 2*j, nx, ny, nxperm, nyperm);
            vy_var(q) = vy_var(q) + 2*(k_/mu/phi)^2 * dpdy * (CYP_p-CYP_n)/(2*del_y);
        end
    end
end
end
