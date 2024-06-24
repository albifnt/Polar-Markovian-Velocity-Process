function vx_mean = calc_vx_mean( Grid, Fluid, Perm, bCond, p0 )
nx = Grid.nx;
ny = Grid.ny;
del_x = Grid.del_x;
mu = Fluid.mu;
phi = Grid.poro;
k_avg = Perm.k_avg;
p0 = reshape(p0,nx,ny)';
vx_mean = zeros((Grid.nx)*Grid.ny,1);
for i=1:nx
    for j=1:ny
        q = i+(j-1)*(nx);
        if( i == 1 )
            if bCond.cond(1) == 1
                k_ = k_avg(j,i);
                vx_mean(q) = 2*k_/mu/phi*(bCond.value(1)-p0(j, i))/del_x;
            else
                vx_mean(q) = bCond.value(1)/phi;
            end
        elseif( i == nx ) 
            if bCond.cond(2) == 1
                k_ = k_avg(j, i);
                vx_mean(q) = 2*k_/mu/phi*(p0(j, i)-bCond.value(2))/del_x;
            else
                vx_mean(q) = -bCond.value(2)/phi;
            end
        else
            k_ = k_avg(j, i); 
            vx_mean(q) = k_/mu/phi*(p0(j,i-1)-p0(j,i+1))/(2*del_x);
        end
    end
end
end
