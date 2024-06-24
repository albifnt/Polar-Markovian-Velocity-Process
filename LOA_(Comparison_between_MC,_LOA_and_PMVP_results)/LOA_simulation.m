% SCRIPT: Small Perturbation Expansion Method
% MASTER STUDENT: Alberto Fontebasso
% PROJECT: Master Thesis
% SUPERVISOR: PD Dr. Meyer-Massetti
% LAB: IFD
clear all

%% COMPUTATIONAL DOMAIN DISCRETIZATION
% Number of cells along the x and y direction of the computational domain
Grid.nx=41;
Grid.ny=41;

% Grid size [m]
Grid.del_x=16/41;
Grid.del_y=16/41;
Grid.del_z=1;

% Number of grid points along the x and y direction for the permeability
% field: to compute the covariances appearing on the right-hand side of SME
% even points at cell interfaces are taken into account 
Grid.nxperm= 2*Grid.nx+1;
Grid.nyperm= 2*Grid.ny+1;

% Grid matrices: these matrices contain the actual discretized location
[Grid.xx, Grid.yy] = meshgrid( linspace(-8 + Grid.del_x/2, 8 - Grid.del_x/2, Grid.nx), linspace(-8 + Grid.del_y/2, 8 - Grid.del_y/2, Grid.ny) );
[Grid.xx_perm, Grid.yy_perm] = meshgrid( linspace(-8, 8, Grid.nxperm), linspace(-8, 8, Grid.nyperm) );



%% GEOSTATISTICAL PROPERTIES OF THE SITE
% path = './Conductivity_16ly_x_16ly_128_x_128/';
Grid.poro = 1;                                          % Porosity

Perm.Y_avg = dlmread('Y_mean.txt');              
Perm.Y_var = dlmread('Y_variance.txt');         

% Perm.Y_avg = dlmread([path 'Y_mean.txt']);              % Mean log-conductivity matrix (Grid.nx, Grid.ny)
% Perm.Y_var = dlmread([path 'Y_variance.txt']);          % Variance lo-conductivity matrix (Grid.nx, Grid.ny)
Perm.k_avg = exp(Perm.Y_avg);                           % Mean conductivity matrix (Grid.nx, Grid.ny)
Perm.k_var = exp(Perm.Y_var);                           % Variance lo-conductivity matrix (Grid.nx, Grid.ny)
Perm.Y_corr = 1;                                        % Correlation length [m]
Perm.CYY = dlmread('C_matrix.txt'); %get_CYYmatrix(Grid, Perm.Y_var, Perm.Y_corr);% Covariance matrix of hydraulic conductivity



%% FLUID PROPERTIES
Fluid.mu=1;                                             % Viscosity



%% BOUNDARY CONDITIONS
% Two vectors are needed for the definition of the boundary conditions: 
% - bCond.cond is a 4-entry vector, whose entries may assume values 0 and 1;
%   when the entry is set to 0, it means that the corresponding boundary
%   condition is a constant flow boundary condition, while if it is set 
%   1, the boundary condition is a connstant pressure one. First entry -->
%   left boundary; Second entry --> right boundary; Third entry --> bottom
%   boundary; Fourth entry --> top boundary
%
% - bCond.value is a 4-entry vector, whose entries assume the value you
%   want to impose for the corresponding boundary condition. First entry -->
%   left boundary; Second entry --> right boundary; Third entry --> bottom
%   boundary; Fourth entry --> top boundary
% 
% REMARK: the constant flow-rate bc assumes positive value if the flow
% comes into the cell
bCond.cond = [0, 0, 1, 1];    
bCond.value = [0, 0, 2, 1];


%% COMPUTATION OF THE MATRIX SHARED BY THE LEFT-HAND SIDE OF ALL EQUATIONS
A = matrix_A( Grid, Fluid, Perm, bCond );
[L, U]=lu(A);                                           % LU decomposition



%% DETERMINATION OF THE FIRST ORDER STATISTICAL MOMENT FOR THE PRESSURE
pressure_RHS = calc_pressure_RHS( Grid, Fluid, Perm, bCond );            
y = L\pressure_RHS;
pressure = U\y;                                         % Mean pressure vector
pressure_ = reshape(pressure, Grid.nx, Grid.ny)';       % Mean pressure matrix



%% DETERMINATION OF THE SECOND ORDER CROSS-COVARIANCE CYP 
tic
CYP = zeros(Grid.nx*Grid.ny, Grid.nxperm*Grid.nyperm);
for i=1:Grid.nxperm
    for j=1:Grid.nyperm
        q = i+(j-1)*Grid.nxperm;
        CYP_RHS = calc_CYP_RHS( Grid, Fluid, Perm, bCond, i, j, pressure);
        y=L\CYP_RHS;
        CYP(:,q)=U\y;
    end
end
toc



%% DETERMINATION OF THE SECOND ORDER STATISTICAL MOMENT CPP FOR THE PRESSURE
tic
CPP = zeros(Grid.nx*Grid.ny, Grid.nx*Grid.ny);
for i=1:Grid.nx
    for j=1:Grid.ny
        q = i+(j-1)*Grid.nx;
        CPP_RHS = calc_CPP_RHS( Grid, Fluid, Perm, bCond, i, j, pressure, CYP);
        y=L\CPP_RHS;
        CPP(:,q)=U\y;
    end
end
toc
for i=1:Grid.nx*Grid.ny
    p_var(i)=CPP(i,i);
end
p_var_=reshape(p_var,Grid.nx,Grid.ny)';



%% DETERMINATION OF THE VELOCITY FIRST ORDER STATISTICAL MOMENTS
% Mean velocity along x direction
vx_mean = calc_vx_mean( Grid, Fluid, Perm, bCond, pressure ); 
vx_mean_ = reshape(vx_mean,Grid.nx,Grid.ny)';

% Mean velocity along y direction
vy_mean = calc_vy_mean( Grid, Fluid, Perm, bCond, pressure ); 
vy_mean_ = reshape(vy_mean,Grid.nx,Grid.ny)';



%% DETERMINATION OF THE VELOCITY SECOND ORDER STATISTICAL MOMENTS
% Variance of velocity along x direction
vx_var = calc_vx_var ( Grid, Fluid, Perm, bCond, pressure, CYP, CPP );
vx_var_ = reshape(vx_var,Grid.nx, Grid.ny)';

% Variance of velocity along y direction
vy_var = calc_vy_var ( Grid, Fluid, Perm, bCond, pressure, CYP, CPP );
vy_var_ = reshape(vy_var,Grid.nx, Grid.ny)';



%% PLOT
% Pressure
figure(1)
set(gcf, 'position', [100 50 1000 500])
subplot(1,2,1)
contourf(Grid.xx, Grid.yy, pressure_, 500, 'LineStyle', 'none');
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Mean pressure');
setCMRcolormap(true);
colorbar;
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.08;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)

subplot(1,2,2)
contourf(Grid.xx, Grid.yy, p_var_, 500, 'LineStyle', 'none');
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Variance pressure');
setCMRcolormap(true);
colorbar;
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.58;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)


% Velocity along the x direction
figure(2)
set(gcf, 'position', [100 50 1000 500])
subplot(1,2,1)
contourf(Grid.xx, Grid.yy, vx_mean_, 500, 'LineStyle', 'none');
axis equal;
axis([-4 4 -4 4])
xlabel('x direction');
ylabel('y direction');
title('Mean velocity x');
setCMRcolormap(true);
colorbar;
caxis([min( min(vx_mean_) ) max( max(vx_mean_) )])
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.08;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)

subplot(1,2,2)
contourf(Grid.xx, Grid.yy, vx_var_, 500, 'LineStyle', 'none');
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Velocity x variance');
setCMRcolormap(true);
colorbar;
caxis([min( min(vx_var_) ) max( max(vx_var_) )])
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.58;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)


% Velocity along the y direction
figure(3)
set(gcf, 'position', [100 50 1000 500])
subplot(1,2,1)
contourf(Grid.xx, Grid.yy, vy_mean_, 500, 'LineStyle', 'none');
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Mean velocity y');
setCMRcolormap(true);
colorbar;
caxis([min( min(vy_mean_) ) max( max(vy_mean_) )])
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.08;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)

subplot(1,2,2)
contourf(Grid.xx, Grid.yy, vy_var_, 500, 'LineStyle', 'none');
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Velocity y variance');
setCMRcolormap(true);
colorbar;
caxis([min( min(vy_var_) ) max( max(vy_var_) )])
set(gca,'FontSize',16)
pos = get(gca, 'Position');
pos(1) = 0.58;
pos(2) = 0.1;
pos(3) = 0.34;
% pos(4) = 0.9;
set(gca, 'Position', pos)



% Mean velocity field
figure(4)
xy = even_stream_arrow(Grid.xx, Grid.yy, vx_mean_, vy_mean_, 3, 3, 'Color', 'k', 'LineWidth',1);
axis equal;
xlabel('x direction');
ylabel('y direction');
title('Mean velocity field');



% %% WRITE TO FILE
% dlmwrite('Mean_x.txt', vx_mean_);
% dlmwrite('Mean_y.txt', vy_mean_);
% dlmwrite('Variance_x.txt', vx_var_);
% dlmwrite('Variance_y.txt', vy_var_);
% dlmwrite('Mean_pressure.txt', pressure_);
% dlmwrite('Variance_pressure.txt', p_var_);



figure(5)
set(gcf, 'position', [100 50 1000 500])
subplot(1,2,1)
contourf(Grid.xx, Grid.yy, Perm.Y_avg, 500, 'LineStyle', 'none')
axis equal
title('Mean Y')
setCMRcolormap(true);
colorbar;


subplot(1,2,2)
contourf(Grid.xx, Grid.yy,Perm.Y_var, 500,'LineStyle', 'none')
axis equal
title('Variance Y')
setCMRcolormap(true);
colorbar;


