% SCRIPT: Velocity variances with interpolating functions
% MASTER STUDENT: Alberto Fontebasso
% PROJECT: Master Thesis
% SUPERVISOR: PD Dr. Meyer-Massetti
% LAB: IFD

clear all
%--------------------------------------------------------------------------%
%% FILE READING MC
% MC
v_x1_mean = dlmread('./velocity_x1_mean.txt');
Mean_y = ( v_x1_mean(1:end-1,:) + v_x1_mean(2:end,:) )*0.5; 
v_x2_mean = dlmread('./velocity_x2_mean.txt');
Mean_x = ( v_x2_mean(:,1:end-1) + v_x2_mean(:,2:end) )*0.5;
[xx, yy] = meshgrid( linspace(-8 + 16/65/2, 8 - 16/65/2, 65), linspace(-8 + 16/65/2, 8 - 16/65/2, 65) );
st_dev_cond_Y = dlmread('./st_dev_cond_Y.txt');


%% ANGLE 
Angle = atan(Mean_x./Mean_y);


%% MEAN VELOCITY FIELD
U = sqrt( Mean_x.^2 + Mean_y.^2 );


%% INTERPOLATED VARIANCES
variance_x = ( 0.14*(st_dev_cond_Y.^2).*exp( 0.289*(st_dev_cond_Y.^2) ) ).* U.^2;
variance_y = ( 0.341*(st_dev_cond_Y.^2).*exp( 0.173*(st_dev_cond_Y.^2) ) ).*U.^2;

VARIANCE_x = zeros(65,65);
VARIANCE_y = zeros(65,65);


%% TRANSFORMATION TO THE DOMAIN COORDINATE SYSTEM
for i=1:length(xx(1,:))
   for j=1:length(yy(:,1))
       COV_STRLN = [variance_x(j,i) 0; 0 variance_y(j,i)]; 
       A = cos(Angle(j,i));
       B = sin(Angle(j,i));
       ROT_MAT = [A -B; B A];
       COV_DOM = ROT_MAT*COV_STRLN*(ROT_MAT');
       VARIANCE_x(j,i) = COV_DOM(1,1);
       VARIANCE_y(j,i) = COV_DOM(2,2);
   end
end


figure(1)
set(gcf, 'position', [100 50 1000 400])
subplot(1,2,1)
contourf(xx, yy, VARIANCE_y, 500, 'LineStyle', 'none')
axis equal
axis([-4 4 -4 4])
setCMRcolormap(true);
colorbar;
caxis([0.5*10^-4 2*10^-4])
xlabel('x_2 / l_Y')
ylabel('x_1 / l_Y')
xticks([-4 -2 0 2 4])
yticks([-4 -2 0 2 4])
set(gca,'FontSize',14)
text('Units', 'normalized', 'Position', [0.05 0.90], 'FontName', 'times', 'String', '(a)','Fontsize', 22)

subplot(1,2,2)
contourf(xx, yy, VARIANCE_x, 500, 'LineStyle', 'none')
axis equal
axis([-4 4 -4 4])
setCMRcolormap(true);
colorbar;
caxis([4*10^-5 8*10^-5])
xlabel('x_2 / l_Y')
ylabel('x_1 / l_Y')
xticks([-4 -2 0 2 4])
yticks([-4 -2 0 2 4])
set(gca,'FontSize',14)
text('Units', 'normalized', 'Position', [0.05 0.90], 'FontName', 'times', 'String', '(b)','Fontsize', 22)

