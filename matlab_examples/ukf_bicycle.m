%% UKF bicycle test
clear all
close all

% load params from file
load('bicycle_data.mat') 

use_laser = 1;
use_radar = 1;

stop_for_sigmavis = false;

%% Data Initialization
x_pred_all = []; % predicted state history
x_est_all = []; % estimated state history with time at row number 6
NIS_radar_all = []; % estimated state history with time at row number 6
NIS_laser_all = []; % estimated state history with time at row number 6

est_pos_error_squared_all = [];
laser_pos_error_squared_all = [];

P_est = 0.2*eye(n_x); % initial uncertainty
P_est(4,4) = 0.3; % initial uncertainty
P_est(5,5) = 0.3; % initial uncertainty

%% process noise

acc_per_sec = 0.3; % acc in m/s^2 per sec
yaw_acc_per_sec = 0.3; % yaw acc in rad/s^2 per sec

Z_l_read = [];

std_las1 = 0.15;
std_las2 = 0.15;

std_radr = 0.3;
std_radphi = 0.03;
std_radrd = 0.3;


% UKF params
n_aug = 7;
kappa = 3-n_aug;

w = zeros(2*n_aug+1,1);
w(1) = kappa/(kappa+n_aug);

for i=2:(2*n_aug+1)
  w(i) = 0.5/(n_aug+kappa);
end

%% UKF filter recursion
%x_est_all(:,1) = GT(:,1);
Xi_pred_all = [];
Xi_aug_all = [];
x_est = [0.1 0.1 0.1 0.1 0.01];
last_time = 0;


% load measurement data from file
fid = fopen('obj_pose-laser-radar-synthetic-ukf-input.txt');

%% State Initialization
tline = fgets(fid); % read first line

% find first laser measurement
while tline(1) ~= 'L' % laser measurement
    tline = fgets(fid); % go to next line
end
    
line_vector = textscan(tline,'%s %f %f %f %f %f %f %f %f %f');
last_time = line_vector{4};
x_est(1) = line_vector{2}; % initialize position p_x
x_est(2) = line_vector{3}; % initialize position p_y

tline = fgets(fid); % go to next line 

% counter 
k = 1;
while ischar(tline)  % go through lines of data file
    
    % find time of measurement
    if tline(1) == 'L' % laser measurement
        if use_laser == false
            tline = fgets(fid); % skip this line and go to next line
            continue;
        else % read laser meas time
            line_vector = textscan(tline,'%s %f %f %f %f %f %f %f %f %f');
            meas_time = line_vector{1,4};
        end
    elseif  tline(1) == 'R' % radar measurement 
        if use_radar == false
            tline = fgets(fid); % skip this line and go to next line
            continue;
        else % read radar meas time
            line_vector = textscan(tline,'%s %f %f %f %f %f %f %f %f %f %f');
            meas_time = line_vector{5};
        end
    else % neither laser nor radar
        disp('Error: not laser nor radar')
        return;
    end
    
    
    delta_t_sec = ( meas_time - last_time ) / 1e6; % us to sec
    last_time = meas_time;

    
    %% Prediction part
    p1 = x_est(1);
    p2 = x_est(2);
    v = x_est(3);
    yaw = x_est(4);
    yaw_dot = x_est(5);
    x = [p1; p2; v; yaw; yaw_dot];
    
    
    std_a = acc_per_sec;     % process noise long. acceleration
    std_ydd = yaw_acc_per_sec;   % process noise yaw acceleration
  
    if std_a == 0;
        std_a = 0.0001;
    end
    if std_ydd == 0;
        std_ydd = 0.0001;
    end
    % Create sigma points
    x_aug = [x ; 0 ; 0];
    P_aug = [P_est zeros(n_x,2) ; zeros(2,n_x) [std_a^2 0 ; 0 std_ydd^2 ]];
    
    %P_aug = nearestSPD(P_aug);
    
    Xi_aug = zeros(n_aug,2*n_aug+1);
    sP_aug = chol(P_aug,'lower');
    Xi_aug(:,1) = x_aug;
    
    for i=1:n_aug
        Xi_aug(:,i+1) = x_aug + sqrt(n_aug+kappa) * sP_aug(:,i);
        Xi_aug(:,i+1+n_aug) = x_aug - sqrt(n_aug+kappa) * sP_aug(:,i);
    end
    
    
    % Predict sigma points
    Xi_pred = zeros(n_x,2*n_aug+1);
    for i=1:(2*n_aug+1)
        p1 = Xi_aug(1,i);
        p2 = Xi_aug(2,i);
        v = Xi_aug(3,i);
        yaw = Xi_aug(4,i);
        yaw_dot = Xi_aug(5,i);
        
        nu_a = Xi_aug(6,i);
        nu_yaw_dd = Xi_aug(7,i);
        
        if abs(yaw_dot) > 0.001
            p1_p = p1 + v/yaw_dot * ( sin (yaw + yaw_dot*delta_t_sec) - sin(yaw));
            p2_p = p2 + v/yaw_dot * ( cos(yaw) - cos(yaw+yaw_dot*delta_t_sec) );
        else
            p1_p = p1 + v*delta_t_sec*cos(yaw);
            p2_p = p2 + v*delta_t_sec*sin(yaw);
        end
        
        v_p = v;
        yaw_p = yaw + yaw_dot*delta_t_sec;
        yaw_dot_p = yaw_dot;
        
        % add noise
        p1_p = p1_p + 0.5*nu_a*delta_t_sec^2 * cos(yaw);
        p2_p = p2_p + 0.5*nu_a*delta_t_sec^2 * sin(yaw);
        v_p = v_p + nu_a*delta_t_sec;
        
        yaw_p = yaw_p + 0.5*nu_yaw_dd*delta_t_sec^2;
        yaw_dot_p = yaw_dot_p + nu_yaw_dd*delta_t_sec;
        
        Xi_pred(1,i) = p1_p;
        Xi_pred(2,i) = p2_p;
        Xi_pred(3,i) = v_p;
        Xi_pred(4,i) = yaw_p;
        Xi_pred(5,i) = yaw_dot_p;
    end
    
    % average and covar of sigma points
    x_pred = 0;
    P_pred = zeros(5,5);
    
    for i=1:2*n_aug+1
        x_pred = x_pred + w(i)* Xi_pred(:,i);
    end
    
    for i=1:2*n_aug+1
        P_pred = P_pred + w(i)* (Xi_pred(:,i) - x_pred)*(Xi_pred(:,i) - x_pred)';
    end
    
    

    %% visualize sigma point examples
    if stop_for_sigmavis && k == 25
        disp('Stopping for sigma point visualization');
        
        % 2d example
        P_s = P_est (1:2,1:2);
        x_s = x(1:2);
        Xi_s = zeros(2,5);
        A = chol(P_s,'lower');
        Xi_s(:,1) = x_s;
        
        for i=1:2
            Xi_s(:,i+1) = x_s + sqrt(3) * A(:,i);
            Xi_s(:,i+1+2) = x_s - sqrt(3) * A(:,i);
        end
        
        error_ellipse(P_s,x_s,'conf', 0.4, 'style', 'k-');
        
        Xi_aug_p1 =  squeeze(Xi_s(1,:,:));
        Xi_aug_p2 =  squeeze(Xi_s(2,:,:));
        
        hold on;
        plot(Xi_aug_p1, Xi_aug_p2, 'or');
        legend('P', 'sigma points')
        axis equal
        
        xlabel('p_x in m');
        ylabel('p_y in m');
        
        save('sigma_visualization.mat', 'x_s','P_s','A','Xi_s', 'Xi_aug', 'Xi_pred');
        
        %return;
    end
    k=k+1;
    
    
    
    %% Update part 
    
    if tline(1) == 'L' % laser measurement
        
        % measurement sigma points
        Xi_z = Xi_pred(1:2,:);
        
        % predict measurement
        z_pred = [0;0];
        for i=1:2*n_aug+1
            z_pred = z_pred + w(i) * Xi_z(:,i);
        end
        
        % measurement covar
        S = zeros(2,2);
        for i=1:2*n_aug+1
            S = S + w(i) * (Xi_z(:,i) - z_pred)* (Xi_z(:,i) - z_pred)';
        end
        S = S + [std_las1^2 0 ; 0 std_las2^2 ];
        
        H = [1 0 0 0 0
            0 1 0 0 0];
        S_lin = H*P_pred*H' + [std_las1^2 0 ; 0 std_las2^2 ];
        
        lindiff = S - S_lin;
        
        % cross correlation
        Tc = zeros(n_x,2);
        for i=1:2*n_aug+1
            Tc = Tc + w(i) * (Xi_pred(:,i) - x_pred) * (Xi_z(:,i) - z_pred)';
        end
        
        % Kalman gain
        K = Tc/S;
        
        K_lin = P_pred* H' / S_lin;
        lindiff2 = K - K_lin;
        
        % update
        z = [line_vector{2} ; line_vector{3} ];
        
        p1_true = line_vector{5};
        p2_true = line_vector{6};
            
        x_est = x_pred + K * (z - z_pred);
        P_est = P_pred - K*S*K'; 
        
        gamma = z - z_pred;
        
        NIS_laser = gamma' * S^(-1) * gamma;
        NIS_laser_all(:,end+1) = [NIS_laser; meas_time];

        % store laser measurement in vector for plotting
        Z_l_read(:,end+1) = [ z ; meas_time];
        
    elseif tline(1) == 'R' % radar measurement
        
        % measurement sigma points
        Xi_z = zeros(3,2*n_aug+1);
        for i=1:2*n_aug+1
            p1 = Xi_pred(1,i);
            p2 = Xi_pred(2,i);
            v  = Xi_pred(3,i);
            yaw = Xi_pred(4,i);
            
            v1 = cos(yaw)*v;
            v2 = sin(yaw)*v;
            
            Xi_z(1,i) = sqrt(p1^2 + p2^2);                       %r
            Xi_z(2,i) = atan2(p2,p1);                             %phi
            Xi_z(3,i) = (p1*v1 + p2*v2 ) / sqrt(p1^2 + p2^2);    %r_dot
        end
        
        % predict measurement
        z_pred = [0;0;0];
        for i=1:2*n_aug+1
            z_pred = z_pred + w(i) * Xi_z(:,i);
        end
        
        % measurement covar
        S = zeros(3,3);
        for i=1:2*n_aug+1
            S = S + w(i) * (Xi_z(:,i) - z_pred)* (Xi_z(:,i) - z_pred)';
        end
        S = S + [std_radr^2  0           0
            0          std_radphi^2  0
            0          0           std_radrd^2];
        
        % cross correlation
        Tc = zeros(n_x,3);
        for i=1:2*n_aug+1
            Tc = Tc + w(i) * (Xi_pred(:,i) - x_pred) * (Xi_z(:,i) - z_pred)';
        end
        
        % Kalman gain
        K = Tc/S;
        
        % update
        z = [line_vector{2};line_vector{3};line_vector{4}];

        p1_true = line_vector{6};
        p2_true = line_vector{7};
        
        gamma = z - z_pred;
        NIS_radar = gamma' * S^(-1) * gamma;
        
        
        
        while (gamma(2)) > pi
            gamma(2) = gamma(2) - 2*pi;
        end
        while (gamma(2)) < -pi
            gamma(2) = gamma(2) + 2*pi;
        end
        
        x_est = x_pred + K * (gamma);
        P_est = P_pred - K*S*K';
        
        NIS_radar_all(:,end+1) = [NIS_radar; meas_time];
        
      
    else
        disp('Error: not laser nor radar')
        return;
    end
    
    err = (x_est(1) - p1_true)^2 + (x_est(2) - p2_true)^2
    
    est_pos_error_squared_all(:,end+1)  =  (x_est(1) - p1_true)^2 + (x_est(2) - p2_true)^2 ;

    x_est_all(:,end+1) = [x_est; meas_time];
    x_pred_all(:,end+1) = x_pred;
    Xi_aug_all(:,:,end+1) = Xi_aug;
    Xi_pred_all(:,:,end+1) = Xi_pred;
    
    
%     figure(3)
%     hold on;
%     plot(GT(1,k), GT(2,k), '-og');
%     plot(x_est(1,:), x_est(2,:), '-or');
%     plot(Z_l(1,k), Z_l(2,k), '-xb');
%     axis equal
%     legend('GT', 'est', 'Lasermeas')
%     k
        
   tline = fgets(fid); % read the next line of the data file
end
fclose(fid);

Xi_pred_p1 =  squeeze(Xi_pred_all(1,:,:));
Xi_pred_p2 =  squeeze(Xi_pred_all(2,:,:));

figure(2)
hold on;
plot(GT(1,:), GT(2,:), '-og'); 
plot(x_est_all(1,:), x_est_all(2,:), '-or');
plot(x_pred_all(1,:), x_pred_all(2,:), '.b');
plot(Xi_pred_p1, Xi_pred_p2, 'xb');
legend('GT', 'est', 'pred', 'Xi pred')

figure(3)
hold on;
plot(GT(1,:), GT(2,:), '-og'); 
plot(x_est_all(1,:), x_est_all(2,:), '-or');
plot(Z_l_read(1,:), Z_l_read(2,:), '-xb');
axis equal
legend('GT', 'est', 'Lasermeas')



%%
figure(1)
hold on;
plot(GT(8,:),GT(1,:), '.-c'); 
plot(x_est_all(6,:),x_est_all(1,:), '-r'); 
plot(Z_l(3,:),Z_l(1,:), '-k'); 
plot(GT(8,:),GT(2,:), '.-b'); 
plot(x_est_all(6,:),x_est_all(2,:), '-r'); 
plot(Z_l(3,:),Z_l(2,:), '-k'); 
plot(GT(8,:),GT(3,:), '.-g');
plot(x_est_all(6,:),x_est_all(3,:), '-g');
plot(GT(8,:),GT(4,:), '.-r'); 
plot(x_est_all(6,:),x_est_all(4,:), '-r'); 
plot(GT(8,:),GT(5,:), '.-m'); 
plot(x_est_all(6,:),x_est_all(5,:), '-m'); 
plot(GT(8,:),[0 diff(GT(3,:))/delta_t_sec], '-c'); 
plot(GT(8,:),[0 diff(GT(5,:))/delta_t_sec], '.c'); 
legend('p1', 'p1est','p1meas', 'p2', 'p2est','p2meas', 'v', 'vest', 'yaw', 'yawest', 'yawrate', 'yawest', 'acc', 'yawacc')

figure(4)
hold on;
if use_radar
plot(NIS_radar_all(2,:) ,NIS_radar_all(1,:)./7.815, '-r'); 
end
if use_laser
plot(NIS_laser_all(2,:),NIS_laser_all(1,:)./5.991, '-b'); 
end
%plot( [1 1], [NIS_laser_all(2,1) NIS_laser_all(2,end)], '-k'); 

legend('NIS radar', 'NIS laser')



RMSE_filter = sqrt(mean(est_pos_error_squared_all))
