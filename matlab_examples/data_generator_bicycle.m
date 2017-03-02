%% UKF bicycle test
clear all
close all

%% Initialization

% set starting position
p_x_start = 0.6;
p_y_start = 0.6;

%TOFILE
delta_t_us = 0.5*1e5; % delta t in microseconds
delta_t_sec = delta_t_us / 1e6; % delta t in sec
round_time_sec = 25; % round time in sec

n_z = round_time_sec/delta_t_sec;  % number of measurements
n_x = 5;  % state dimension

%TOFILE
GT = zeros(n_x,n_z); % ground truth: p1 p2 v_abs yaw yaw_dot v1 v2

% set start values
GT(1,1) =  p_x_start;
GT(2,1) =  p_y_start;

%TOFILE
Z_l = zeros(3,n_z); % laser measurements: pos1, pos2 time                 
%TOFILE
Z_r = zeros(4,n_z); % radar measurements: r, phi, r_dot time

std_las1 = 0.15;
std_las2 = 0.15;

std_radr = 0.3;
std_radphi = 0.03;
std_radrd = 0.3;

% generate motion example velocity and yaw rate
base_timestamp = 1477010443000000; % October 2016

for k = 1:n_z
    timestamp = base_timestamp+(k-1)*delta_t_us;
    
    % time scaling
    round_time_sec = 25; % round time in sec
        
    % timestamps
    GT(8,k) = timestamp;
    % velocity
    time_sec = k*delta_t_us/1e6;
    GT(3,k) = 2 + 0.2*cos(2*2*pi/round_time_sec*time_sec);
    % yaw rate
    GT(5,k) = 0.55*sin(2*pi/round_time_sec*time_sec);
end

Z_l(3,1)=base_timestamp;
Z_r(4,1)=base_timestamp;


% first laser measurement (special case) 
Z_l(1,1) = GT(1,1) + normrnd(0, std_las1) ;
Z_l(2,1) = GT(2,1) + normrnd(0, std_las2) ;
Z_l(3,1) = GT(8,1) ;  % time

% radar measurements
p1 = GT(1,1);
p2 = GT(2,1);
v = GT(3,1);
yaw = GT(4,1);

v1 = GT(3,1) * cos(GT(4,1));
v2 = GT(3,1) * sin(GT(4,1));

Z_r(1,1) = sqrt(p1^2 + p2^2) + normrnd(0, std_radr) ;  %r
Z_r(2,1) = atan2(p2,p1) + normrnd(0, std_radphi) ;                           %phi
Z_r(3,1) = (p1*v1 + p2*v2 ) / sqrt(p1^2 + p2^2) + normrnd(0, std_radrd) ;  %r_dot
Z_r(4,1) = GT(8,1) ;  % time

GT(6,1) = cos(yaw)* v; % v_x
GT(7,1) = sin(yaw)* v; % v_y


% positions and measurements
for k = 2:n_z
    
    delta_t_s = (GT(8,k) - GT(8,k-1))/1e6;
    
    p1 = GT(1,k-1);
    p2 = GT(2,k-1);
    v = GT(3,k-1);
    yaw = GT(4,k-1);
    yaw_dot = GT(5,k-1);
    
    if abs(yaw_dot) > 0.001   
        p1_p = p1 + v/yaw_dot * ( sin (yaw + yaw_dot*delta_t_s) - sin(yaw)); 
        p2_p = p2 + v/yaw_dot * ( cos(yaw) - cos(yaw+yaw_dot*delta_t_s) );
    else
        p1_p = p1 + v*delta_t_s*cos(yaw);
        p2_p = p2 + v*delta_t_s*sin(yaw);
    end
    v_p = v;
    yaw_p = yaw + yaw_dot*delta_t_s;
    yaw_dot_p = yaw_dot;
    
    GT(1,k) = p1_p;
    GT(2,k) = p2_p;
    GT(4,k) = yaw_p;

    % laser measurements    
    Z_l(1,k) = GT(1,k) + normrnd(0, std_las1) ;
    Z_l(2,k) = GT(2,k) + normrnd(0, std_las2) ;
    Z_l(3,k) = GT(8,k) ;  % time

    % radar measurements
    p1 = GT(1,k);
    p2 = GT(2,k);
    v = GT(3,k);
    yaw = GT(4,k);

    v1 = GT(3,k) * cos(GT(4,k));
    v2 = GT(3,k) * sin(GT(4,k));

    Z_r(1,k) = sqrt(p1^2 + p2^2) + normrnd(0, std_radr) ;  %r
    Z_r(2,k) = atan2(p2,p1) + normrnd(0, std_radphi) ;                           %phi
    Z_r(3,k) = (p1*v1 + p2*v2 ) / sqrt(p1^2 + p2^2) + normrnd(0, std_radrd) ;  %r_dot
    Z_r(4,k) = GT(8,k) ;  % time
    
    GT(6,k) = cos(yaw)* v; % v_x
    GT(7,k) = sin(yaw)* v; % v_y
end

%% GENERATE UKF GROUND 
GTT=GT';

%GT - p1 p2 v_abs yaw yaw_dot v1 v2
fileID = fopen('obj_pose-laser-radar-synthetic-gt.txt','w');
% nbytes = fprintf(fileID,'%d\t%d\n',GTT);
% fclose(fileID);

[nrows,ncols] = size(GTT);
for row = 1:nrows
    fprintf(fileID,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',GTT(row,1:8));
end
fclose(fileID);

%% WRITE EKF/UKF INPUT MEASUREMENT FILE
%add timestamps
%L    p1    p2    timestamp=
%R    r theta    r_dot    timestamp

%laser structure
Z_lT = Z_l';
N_L =size(Z_lT,1);
v_z = zeros(N_L,1); %vector of zeros
Z_lT = [Z_lT, v_z];
%Z_lT = [Z_lT];
%create a laser cell array - mixed numbers and characters
C_l = cell(N_L,1);
C_l(:) = {'L'};
C_l = [ C_l, num2cell(Z_lT)];

%radar structure
Z_rT = Z_r';
N_R =size(Z_rT,1);
%create a radar cell array - mixed numbers and characters
C_r = cell(N_R,1);
C_r(:) = {'R'};
C_r = [ C_r, num2cell(Z_rT)];


AB ={};
for i = 1:size(Z_lT,1)
    %append elements from both matrices
    AB = [AB; C_l(i,:)];
    AB = [AB; C_r(i,:)];
  
end



fileID = fopen('obj_pose-laser-radar-synthetic-ukf-input.txt','w');
[nrows,ncols] = size(AB);

row_gt = int16(1);
rm_index = 0;
for row = 1:nrows
    %ismember( 'L', AB{row,1} );
    rm_index = rm_index + 1;
    rm_index = mod(rm_index, 4)
    if( AB{row,1} == 'L' &&  mod(row,4) == 1 )
        
        %when 
        %estimations
        fprintf(fileID,'%s\t%d\t%d\t%d\t',AB{row,1:ncols-1});
        %ground truth speed and position
        fprintf(fileID,'%d\t%d\t%d\t%d\t',GTT(row_gt, [1 2 6 7]));
        
        %yaw and yaw_rate
         fprintf(fileID,'%d\t%d\n',GTT(row_gt, [4 5]));
        
    elseif( AB{row,1} == 'R' &&  mod(row,4) == 0 )
        
        %estimations
        fprintf(fileID,'%s\t%d\t%d\t%d\t%d\t',AB{row,:});
        %ground truth speed and position
        fprintf(fileID,'%d\t%d\t%d\t%d\t',GTT(row_gt, [1 2 6 7]));
         %yaw and yaw_rate
         fprintf(fileID,'%d\t%d\n',GTT(row_gt, [4 5]));
         
    end
    if(mod(row,2) == 0)
        row_gt = row_gt + 1;
    end;
end

% nbytes = fprintf(fileID,'%s\t%d\t%d\t%d\t%d\t\n',AB);
fclose(fileID); 

% save the workspace in file
save('bicycle_data.mat') 

figure(2)
hold on;
plot(GT(1,:), GT(2,:), '-og'); 
xlabel('x');
ylabel('y');
axis equal
legend('GT')

%%
figure(1)
hold on;
plot(GT(2,:), '.-b'); 
plot(GT(3,:), '.-g');
plot(GT(4,:), '.-r'); 
plot(GT(5,:), '.-m'); 
plot(diff(GT(3,:))/delta_t_sec, '-c'); 
plot(diff(GT(5,:))/delta_t_sec, '.c'); 
legend('p2', 'v', 'yaw', 'yawrate', 'acc', 'yawacc')


% TODO: fill file with consecutive time stamps, fill GT into file.
