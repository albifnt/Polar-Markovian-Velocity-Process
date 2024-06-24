% SCRIPT: MC simulation
% MASTER STUDENT: Alberto Fontebasso
% PROJECT: Master Thesis
% SUPERVISOR: PD Dr. Meyer-Massetti
% LAB: IFD

clear all
%% READING FILE
files = dir('./Realization_0/*.txt');
count=1;
filename = {};
for i=1:length(files)
    if ~any(strcmp(filename, files(i).name(8:10)))
        filename{count} = files(i).name(8:10);
        count = count + 1;
    end
end
clear files count i

%% GRID SIZE
delta = 2.85;
n = 0.31;

%% ELABORATE DATA
for q=1:18
    
    matrix_U1 = zeros(15,256,100000);
    matrix_U2 = zeros(15,256,100000);
    count = 1;
    for k=0:99
        for i=1:length(filename)

            Result = dlmread(['./Realization_' num2str(k) '/Result_' filename{i} '.txt']);

            % GRID PARAMETERS FOR THE HORIZONTAL
            nx = 272;
            ny = 256;

            % MC RESULTS
            qx = Result(3:((nx+1)*ny+2),1);
            qx = reshape(qx,nx+1,ny)';

            qy = Result(((nx+1)*ny+4):((nx+1)*ny+3+nx*(ny+1)),1);
            qy = reshape(qy,nx,ny+1)';

            h = Result(((nx+1)*ny+3+nx*(ny+1)+2):end,1);
            h = reshape(h,nx,ny)';


            % GRID PARAMETERS FOR THE VERTICAL
            nx = 256;
            ny = 272;

            % MC RESULTS
            qx1 = qx';
            qx2 = qy';
            h = h';

            % MC RESULTS
            ux1 = qx1/(delta*n);
            ux2 = qx2/(delta*n);

            % PLOT
            U_x2 = ( ux2(:,1:end-1) + ux2(:,2:end) )*0.5;
            U_x1 = ( ux1(1:end-1,:) + ux1(2:end,:) )*0.5;

            matrix_U1(:,:,count) = U_x1( (15*q - 14):(15*q), : );
            matrix_U2(:,:,count) = U_x2( (15*q - 14):(15*q), : );
            count = count + 1;
            count
        end
    end

    Media_U1 = zeros(15,256);
    Media_U2 = zeros(15,256);
    Variance_U1 = zeros(15,256);
    Variance_U2 = zeros(15,256);

    for w=1:length(matrix_U1(:,1,1))
       for z=1:length(matrix_U1(1,:,1))

           Media_U1(w,z) = mean(matrix_U1(w,z,:));
           Media_U2(w,z) = mean(matrix_U2(w,z,:));
           Variance_U1(w,z) = var(matrix_U1(w,z,:));
           Variance_U2(w,z) = var(matrix_U2(w,z,:));

       end
    end
    
    dlmwrite(['./MATRICES/Media_U1_' num2str(q) '.txt'], Media_U1);
    dlmwrite(['./MATRICES/Media_U2_' num2str(q) '.txt'], Media_U2);
    dlmwrite(['./MATRICES/Variance_U1_' num2str(q) '.txt'], Variance_U1);
    dlmwrite(['./MATRICES/Variance_U2_' num2str(q) '.txt'], Variance_U2);

end

