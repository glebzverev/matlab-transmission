clear all
close all

data_t15 = load("T15MD_imp_898_1800_135_05_al110-4.dat");
data_var1 = load("poloidal_Flux.txt");
r_t15 = data_t15(1,2:end);
r_var1 = data_var1(1,2:end);
z_t15 = data_t15(2:end,1);
z_var1 = data_var1(2:end,1);

data_t15 = data_t15(2:end,2:end);
data_var1 = data_var1(2:end, 2:end);
delta_r_files = sum(abs(r_t15 - r_var1));
delta_z_files = sum(abs(z_t15 - z_var1));
F = data_t15;
F2 = data_var1;
R = r_t15; Z = z_t15;
dr = R(2)-R(1);
dz = Z(2)-Z(1);
%%
camera_cell = readcell('камера.xlsx','Sheet','camera');
r_camera = [camera_cell{4:13,2}];
z_camera = [camera_cell{4:13,3}];
xl = xlsread('камера.xlsx','camera');
pf_coils = xlsread('камера.xlsx','PF-coils');
pf_coils_cells = readcell('камера.xlsx','Sheet','PF-coils');
inductor = xlsread('камера.xlsx','Inductor');%только числа читает
r_inductor = [inductor(:,2);inductor(:,6);inductor(:,10)];
z_inductor = [inductor(:,3);inductor(:,7);inductor(:,11)];
probes_cell = readcell('zondu.xlsx');%%зонды меняю
coils_cell = readcell('камера.xlsx','Sheet','coils');
probes = zeros(37, 2);%
coils = zeros(27-17+1, 4-3+1);

for i=1:37%
    for j=1:2%
        probes(i,j) = str2num(probes_cell{i,j});
    end
end
for i=17:27
    for j=3:4
        coils(i-16,j-2) = str2num(coils_cell{i,j});
    end
end

r_pf_coils = pf_coils(:,1);
z_pf_coils = pf_coils(:,2);
% dr_pf_coils = pf_coils(:,4);
% dz_pf_coils = pf_coils(:,5);
dr_pf_coils = pf_coils(:,4);
dz_pf_coils = pf_coils(:,5);


r_div_plast = [camera_cell{4:7,18}]';
z_div_plast = [camera_cell{4:7,19}]';
d_div_plast = [camera_cell{4:7,20}]';


r_camera_inner = [camera_cell{16:25,2}]';
z_camera_inner = [camera_cell{16:25,3}]';
r_camera_outer = [camera_cell{16:25,4}]';
z_camera_outer = [camera_cell{16:25,5}]';

[RR, ZZ] = meshgrid(R,Z);%из одномерных делает двумерные
%%
sensors_cell = readcell('sensors.xlsx');%сенсоры
sensors = zeros(30, 2);
for i=1:30%
    for j=1:2%
        sensors(i,j) = str2num(sensors_cell{i,j});
    end
end
%%
coils1_cell = readcell('coils.xlsx');%coils
coils1 = zeros(10, 2);
for i=1:10%
    for j=1:2%
        coils1(i,j) = str2num(coils1_cell{i,j});
    end
end

%%
diafragma = readcell('difragma.xlsx');%сенсоры
diaf = zeros(20, 2);
for i=1:20%
    for j=1:2%
        diaf(i,j) = str2num(diafragma{i,j});
    end
end
%
%%
figure(2)
clf
hold on
grid on
daspect([1 1 1])%маштаб

pgon_camera = polyshape({r_camera_inner,r_camera_outer},...
    {z_camera_inner, z_camera_outer});


plot(pgon_camera, 'FaceColor','m')%рисуем




for i=1:3
    if i == 2
           rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','b') 
      text(pf_coils_cells{2+i,2}+0.1,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
    else
   rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','b') 
      text(pf_coils_cells{2+i,2}+0.1,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
    end
end


for i=4:4
   rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','g') 
      text(pf_coils_cells{2+i,2}+0.1,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
end

for i=5:length(r_pf_coils)-1
   rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','b') 
   
   text(pf_coils_cells{2+i,2}+0.1,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
    
end


for i=6:length(r_pf_coils)-2
   rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','b') 
   
   text(pf_coils_cells{2+i,2}+0.1,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
    
end


for i=length(r_pf_coils):length(r_pf_coils)
   rectangle('Position',[r_pf_coils(i)-dr_pf_coils(i)/2,...
       z_pf_coils(i)-dz_pf_coils(i)/2,dr_pf_coils(i),dz_pf_coils(i)],...
       'FaceColor','g') 
   
   text(pf_coils_cells{2+i,2}+0.26,pf_coils_cells{2+i,3},pf_coils_cells{2+i,1},...
       'FontSize',8, 'HorizontalAlignment','left')
    
end

for i=1:length(probes)
    r_loc = probes(i,1);
    z_loc = probes(i,2);
    rad = 0.02;
    rectangle('Position', [r_loc-rad, z_loc-rad, 2*rad, 2*rad], 'FaceColor', 'black');
end


for i=1:length(sensors)
    r_loc = sensors(i,1);
    z_loc = sensors(i,2);
    rad = 0.02;
    rectangle('Position', [r_loc-rad, z_loc-rad, 2*rad, 2*rad], 'FaceColor', 'yellow');
end

r_d = zeros(1,20);
z_d = zeros(1, 20);
for i=1:length(diaf)
    r_d(i) = diaf(i,1);
    z_d(i) = diaf(i,2);
    %rad = 0.02;
   % rectangle('Position', [r_loc-rad, z_loc-rad, 2*rad, 2*rad], 'FaceColor', 'red');
end

plot(r_d, z_d, 'black', 'LineWidth',1.5)


%for i=1:length(coils1)
 %   r_loc = coils1(i,1);
 %   z_loc =  coils1(i,2);
  %  rad = 0.02;
  %  rectangle('Position', [r_loc-rad, z_loc-rad, 2*rad, 2*rad],
  %  'FaceColor', 'red');%
%end


% for i=1:length(coils)
%     r_loc = coils(i,1);
%     z_loc = coils(i,2);
%     rad = 0.05;
%     rectangle('Position', [r_loc-rad, z_loc-rad, 2*rad, 2*rad], 'Curvature', [1,1], 'FaceColor', 'g');
% end

R=0.776589; 
Z1=0.880668;
Z2=1.08082;
t = 0.01;
%R_11=0.776589; Z_1 = -0.848644; Z_2=-1.0568;
% 3. R1=1.8259, Z1=-1.20091, R2=1.93747, Z2=-1.03278
% 4. R1=1.93747, Z1=1.04079, R2=1.8259, Z2=1.19291


rectangle('Position', [0.776589-t, 0.880668, 2*t, 1.08082 - 0.880668], 'FaceColor', 'c');
rectangle('Position', [0.776589-t, -1.0568, 2*t, abs(-0.848644 - -1.03278)], 'FaceColor', 'c');
t = 0.02;
pgon = polyshape([1.8259-t+0.01 1.93747-t 1.93747+t-0.01 1.8259+t],[-1.20091+0.01 -1.03278 -1.03278-0.01 -1.20091]);
plot(pgon, 'FaceColor', 'c')
pgon = polyshape([1.93747-t 1.8259-t+0.01 1.8259+t 1.93747+t-0.01],[1.04079 1.19291-0.01 1.19291 1.04079+0.01]);
plot(pgon, 'FaceColor', 'c')

%rectangle('Position', [1.8259-t, -1.20091, 2*t, abs(-0.848644 - -1.03278)], 'FaceColor', 'b');
%rectangle('Position', [0.776589-t, -1.0568, 2*t, abs(-0.848644 - -1.03278)], 'FaceColor', 'b');

rectangle('Position', [r_div_plast(1), z_div_plast(1)-d_div_plast(1)/2, r_div_plast(2)-r_div_plast(1), d_div_plast(1)], 'FaceColor', [.7 .7 .7]);
rectangle('Position', [r_div_plast(3), z_div_plast(3)-d_div_plast(3)/2, r_div_plast(4)-r_div_plast(3), d_div_plast(3)], 'FaceColor', [.7 .7 .7]);



%[i_12, j_12, A_min_grad, i_22, j_22, A_minmax] = indices_grad_and_minimax(F2,dr, dz);
%sep = A_minmax(1);
sep = 0.0285771;
%plot(R(i_22),Z(j_22), '*r')
x = -0.53256
;
%[~, hCont0.005755] = contour(RR,ZZ,F2, [sep,sep],'r');

levels = linspace(x,13.5*x,75);
[~, hCont7] = contour(RR,ZZ,F, levels);

x = -0.53256
;
[~, hCont4] = contour(RR,ZZ,F, [x, x],'g');

contourLegend([hCont4], {strcat(num2str(x),' inverse')})


%0.2706E-01

hold off
