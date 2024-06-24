function vx_var = calc_vx_var( Grid, Fluid, Perm, bCond, p0, CYP, CPP )
nx = Grid.nx;
ny = Grid.ny;
nxperm = Grid.nxperm;
nyperm = Grid.nyperm;
del_x = Grid.del_x;
phi = Grid.poro;
mu = Fluid.mu;
k_avg = Perm.k_avg;
p0 = reshape(p0,nx,ny)';
vx_var = zeros((nx)*ny,1);
CYY = Perm.CYY;
for i=1:nx
    for j=1:ny
        q = i+(j-1)*(nx);
        if( i == 1 )
            k_ = k_avg(j, i);
            if bCond.cond(1) == 1 
                dpdx = 2*(p0(j,i) - bCond.value(1))/(del_x);
                vx_var(q) = vx_var(q) + (k_/mu/phi)^2 * dpdx^2 * get_CYY(CYY, 2*i-1, 2*j, 2*i-1, 2*j, nxperm, nyperm);
            else
                vx_var(q)=0;
            end
        elseif( i == nx ) 
            k_ = k_avg(j, i);
            if bCond.cond(2) == 1
                dpdx = 2*(bCond.value(2)-p0(j,i))/(del_x);
                vx_var(q) = vx_var(q) + (k_/mu/phi)^2 * dpdx^2 * get_CYY(CYY, 2*i+1, 2*j, 2*i+1, 2*j, nxperm, nyperm);
            else
                vx_var(q) = 0;
            end
        else
            k_ = k_avg(j, i);
            qz = i + (j-1) * nx; 
            qn = (i-1) + (j-1) * nx;
            qp = i+1 + (j-1) * nx;
            dcpp_n = (CPP(qz,qn) - CPP(qn,qn))/del_x;
            dcpp_p = (CPP(qp,qp) - CPP(qz,qp))/del_x;
            ddcpp = (dcpp_p - dcpp_n)/(2*del_x);
            vx_var(q) = vx_var(q) + (k_/mu/phi)^2 * ddcpp;
            
            dpdx = (p0(j,i+1) - p0(j,i-1))/(2*del_x); %why the derivative at the cell center?
            
            
            vx_var(q) = vx_var(q) + (k_/mu/phi)^2 * dpdx^2 * get_CYY(CYY,2*i,2*j,2*i,2*j, nxperm, nyperm);
            CYP_p = get_CYP(CYP, i+1, j, 2*i, 2*j, nx, ny, nxperm, nyperm);
            CYP_n = get_CYP(CYP, i-1, j, 2*i, 2*j, nx, ny, nxperm, nyperm);
            vx_var(q) = vx_var(q) + 2*(k_/mu/phi)^2 * dpdx * (CYP_p-CYP_n)/(2*del_x);
        end
    end
end
end
