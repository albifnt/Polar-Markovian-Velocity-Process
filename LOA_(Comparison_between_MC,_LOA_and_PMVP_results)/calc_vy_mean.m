function vy_mean = calc_vy_mean( Grid, Fluid, Perm, bCond, p0 )
nx = Grid.nx;
ny = Grid.ny;
del_y = Grid.del_y;
mu = Fluid.mu;
phi = Grid.poro;
k_avg = Perm.k_avg;
p0 = reshape(p0,nx,ny)';
vy_mean = zeros((Grid.nx)*Grid.ny,1);
for i=1:nx
    for j=1:ny
        q = i+(j-1)*(nx);
        if( j == 1 )
            if bCond.cond(3) == 1
                k_ = k_avg(j,i);
                vy_mean(q) = 2*k_/mu/phi*(bCond.value(3)-p0(j, i))/del_y;
            else
                vy_mean(q) = bCond.value(3)/phi;
            end
        elseif( j == ny ) 
            if bCond.cond(4) == 1
                k_ = k_avg(j, i);
                vy_mean(q) = 2*k_/mu/phi*(p0(j, i)-bCond.value(4))/del_y;
            else
                vy_mean(q) = -bCond.value(4)/phi;
            end
        else
            k_ = k_avg(j, i); 
            vy_mean(q) = k_/mu/phi*(p0(j-1,i)-p0(j+1,i))/(2*del_y);
        end
    end
end
end
