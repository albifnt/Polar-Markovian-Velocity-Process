% SCRIPT: Covariance Matrix Evaluation
% MASTER STUDENT: Alberto Fontebasso
% PROJECT: Master Thesis
% SUPERVISOR: PD Dr. Meyer-Massetti
% LAB: IFD
clear all;

nx = 83;
ny = 83;
x = linspace(-8, 8, nx);
y = linspace(-8, 8, ny);
[X,Y] = meshgrid(x,y);
l_Y = 1;
Variance = dlmread('Y_variance_matrix.txt');
C = zeros(nx*ny, nx*ny);

for i=1:nx*ny
    i_i = fix(i/nx) + 1;
    if ( mod(i,nx) == 0)
        i_i = fix(i/nx); 
    end
    j_i = mod(i,nx);
    if ( j_i == 0)
       j_i = nx; 
    end
    for j=i:nx*ny
        i_j = fix(j/nx) + 1;
        if (  mod(j,nx) == 0)
            i_j = fix(j/nx); 
        end
        j_j = mod(j,nx);
        if ( j_j == 0)
            j_j = nx; 
        end
        C(i,j) = sqrt( Variance(i_i, j_i) )*sqrt( Variance(i_j, j_j) )*exp( -sqrt( ( ( X(1, j_i) - X(1, j_j) )/l_Y )^2 + ( ( Y(i_i, 1) - Y(i_j, 1) )/l_Y )^2 ) );
    end
end

C_lower = triu(C,1)';
C = C + C_lower;
dlmwrite('C_matrix.txt', C);