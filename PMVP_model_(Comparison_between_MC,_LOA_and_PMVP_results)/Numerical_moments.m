% SCRIPT: PMVP statistical moments of velocity
% MASTER STUDENT: Alberto Fontebasso
% PROJECT: Master Thesis
% SUPERVISOR: PD Dr. Meyer-Massetti
% LAB: IFD

clear all
%% READING FILE
path = './PMVP_model_(PMVP_particle_generation)/FILE/';
files = dir([path '*.txt']);
count=1;
filename = {};
for i=1:length(files)
    if ~any(strcmp(filename, files(i).name(1)))
        filename{count} = files(i).name(1);
        count = count + 1;
    end
end

%% NUMBER NODES PER DIMENSION = 9
n = 17;
transverse = linspace(-4,4,n);
longitudinal = linspace(-4,4,n);
Count_matrix = zeros(n,n);
Sum_matrix_x1 = zeros(n,n);
Sum_matrix_x2 = zeros(n,n);
Velocity = {};
for i=1:n
   for k=1:n
       Velocity{i,k}.x1 = [];
       Velocity{i,k}.x2 = [];
   end
end

for i=1:length(filename)
    x1 = dlmread([path filename{i} '_x1.txt']);
    x2 = dlmread([path filename{i} '_x2.txt']);
    velocity_x1 = dlmread([path filename{i} '_velocity_x1.txt']);
    velocity_x2 = dlmread([path filename{i} '_velocity_x2.txt']);

    for id_trial = 1:length(x1(:,1))
        for id_time = 1:length(x1(1,:))-1
            
            if ( x1(id_trial,id_time) > 4.5)
                break
            end
            
            [d_x1,id_x1] = min( abs( longitudinal - x1(id_trial,id_time) ) );
            [d_x2,id_x2] = min( abs( transverse - x2(id_trial,id_time) ) );
            if (d_x1 <= 0.5 & d_x2 <= 0.5)
                Count_matrix(id_x1, id_x2) = Count_matrix(id_x1, id_x2) + 1;
                Sum_matrix_x1(id_x1, id_x2) = Sum_matrix_x1(id_x1, id_x2) + velocity_x1(id_trial, id_time + 1);
                Sum_matrix_x2(id_x1, id_x2) = Sum_matrix_x2(id_x1, id_x2) + velocity_x2(id_trial, id_time + 1);
                Velocity{id_x1, id_x2}.x1(end+1) = velocity_x1(id_trial, id_time + 1);
                Velocity{id_x1, id_x2}.x2(end+1) = velocity_x2(id_trial, id_time + 1);
            end
            
        end
        
    end
end

[xx, yy] = meshgrid(transverse, longitudinal);
Velocity_matrix_x1 = Sum_matrix_x1./Count_matrix;
Velocity_matrix_x2 = Sum_matrix_x2./Count_matrix;

Velocity_variance_x1 = zeros(n,n);
Velocity_variance_x2 = zeros(n,n);

for i=1:n
   for k=1:n
       Velocity_variance_x1(i,k) = var(Velocity{i,k}.x1);
       Velocity_variance_x2(i,k) = var(Velocity{i,k}.x2);
   end
end

figure(1)
contourf(xx,yy, Velocity_matrix_x1, 500, 'LineStyle', 'none')
axis equal

figure(2)
contourf(xx,yy, Velocity_matrix_x2, 500, 'LineStyle', 'none')
axis equal

figure(3)
contourf(xx,yy, Velocity_variance_x1, 500, 'LineStyle', 'none')
axis equal

figure(4)
contourf(xx,yy, Velocity_variance_x2, 500, 'LineStyle', 'none')
axis equal

% print(3, '-dpng', '3.png', '-r600');
% print(4, '-dpng', '4.png', '-r600');
% 
dlmwrite('XX.txt',xx);
dlmwrite('YY.txt',yy);
dlmwrite('Mean_U1.txt',Velocity_matrix_x1);
dlmwrite('Mean_U2.txt',Velocity_matrix_x2);
dlmwrite('Variance_U1.txt',Velocity_variance_x1);
dlmwrite('Variance_U2.txt',Velocity_variance_x2);
