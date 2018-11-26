close all
clear all
%  INPUT. The file requires an k x n matrix, i.e., where k is the number of 
%  joints or cartesian coordinates and n are all the points (initial, 
%  intermidiate and final), an vector of size n-1 that specifies the 
%  duration of each segment, the increment of time within each segment,
%  e.g., step=0.01s, and the type of motion ('cyclic' or 'prescribed' [V0,
%  Vf]).

% joints = [1     2   3   4   5   6
%           7     8   9   10  11  12
%           13    14  15  16  17  18
%           19    20  21  22  23  24];

P_ee = [0 500 0 0;                %   x
        0 500 0 0;                %   y
        2000 2000 2000 2000;      %   z
        0 45 45 0;                %   a 
        0 0 0 0;                  %   b   
        0 0.7068766 0.6332458 0]; %   y

tf = [0.40385, 0.8077, 0.8077];


%%  Inverse Kinematic Fixed Variables
%Characteristics of the manipulator
alpha = [-90,30,150]*pi/180;  % Definition of base platform
L = 1500;                     % [mm] base length (centre to joint)
LL = 1000;                    % mobile length (centre to joint)
P = 3;                        % Link first length
PP = 2.5;                     % Length second Length

% Base Angles
Beta = [225
         315
         345
         75
         105
         195];

% Beta = angle to get to X axis for zero displacement config      
Betad = Beta + 90;

%Platform Angles
anglePlat = [255
             285
             15
             45
             135
             165];
         

p0pi = LL*[cosd(anglePlat(1)),cosd(anglePlat(2)),cosd(anglePlat(3)),cosd(anglePlat(4)),cosd(anglePlat(5)),cosd(anglePlat(6));...
          sind(anglePlat(1)),sind(anglePlat(2)),sind(anglePlat(3)),sind(anglePlat(4)),sind(anglePlat(5)),sind(anglePlat(6));
          0,0,0,0,0,0];
             
%%  Inverse Kinematics Loop
for j =1:length(P_ee(1,:))
    % Solve the INVERSE KINEMATICS for the jth posture and store
    % the joint displacements in a (18x1)vector, say joints_j, which
    % includes the joint displacements for branches 1 through 6:
    % joints_i =[theta1_1 theta2_1, L_1, theta1_2 theta2_2, L_2,... ,
    % theta1_6 theta2_6, L_6];
    % Store vectors in matrix form for each end-effector’s position
    
    x0 = P_ee(1,j);
    y0 = P_ee(2,j);
    z0 = P_ee(3,j);
    alpha = P_ee(4,j); 
    beta = P_ee(5,j);
    gamma = P_ee(6,j);
    
    joints_j = [];
    
    %   Yaw, Pitch, Roll Matrix
    R0_ee = [cosd(alpha)*cosd(beta), cosd(alpha)*sind(beta)*sind(gamma)-sind(alpha)*cosd(gamma), cosd(alpha)*sind(beta)*cosd(gamma)+sind(alpha)*sind(gamma)
             sind(alpha)*cosd(beta), sind(alpha)*sind(beta)*sind(gamma)+cosd(alpha)*cosd(gamma), sind(alpha)*sind(beta)*cosd(gamma)-cosd(alpha)*sind(gamma)
             -sind(beta)           , cosd(beta)*sind(gamma)                                    , cosd(beta)*cosd(gamma)];

    for i = 1:6
        %   Generate the position vector of spherical joint for IK
        P = [x0; y0; z0] + LL*R0_ee*[cosd(anglePlat(i));sind(anglePlat(i));0];

        %   Theta 1 Inverse Kinematics equations
        theta1 = atan2d(P(3),P(1)*cosd(Betad(i))+P(2)*sind(Betad(i)));

        %   New prismatic joint length from IK 
        Lnew = sqrt((P(1)*sind(Betad(i))-P(2)*cosd(Betad(i))-L)^2 + (P(3)*sind(theta1) + P(1)*cosd(Betad(i))*cosd(theta1) + P(2)*sind(Betad(i))*cosd(theta1))^2);
       
        %   Theta 2 Inverse Kinematic Equations used as input for the FK
        theta2 = atan2d((P(3)*sind(theta1)+P(1)*cosd(Betad(i))*cosd(theta1)+P(2)*sind(Betad(i))*cosd(theta1))/Lnew,(P(1)*sind(Betad(i))-P(2)*cosd(Betad(i))-L)/Lnew);
        joints_j = [joints_j, theta1, theta2, Lnew]; 
    end
    joints(:,j) = joints_j;
end

%% TRAJECTORY GENERATION. Trajectory generation of all the joints
dt = 0.1; %stepsize
[position,velocity,acceleration,time] = via_points_match_VA(joints, tf, dt, 'prescribed',[0,0]);

%% Plotting Graphs
[n, j] = size(acceleration);    % rows, columns

for i = 1:(n/3)
   %%   Prismatic Plots
   % Plotting Displacement Graph
   hold on
   figure(1)
   plot(time(3*i,:), position(3*i,:)./10000)
   title('Prismatic - Time vs Displacement','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Prismatic Joint Displacement (m)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Velocity Graph
   hold on
   figure(2)
   plot(time(3*i,:), velocity(3*i,:)./10000)
   title('Prismatic - Time vs Velocity','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Prismatic Joint Velocity (m/s)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Acceleration Graph
   hold on
   figure(3)
   plot(time(3*i,:), acceleration(3*i,:)./10000)
   title('Prismatic - Time vs Acceleration','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Prismatic Acceleration (m/s^2)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   %%   Theta 1 Plots
   % Plotting Displacement Graph
   hold on
   figure(4)
   plot(time(3*i-2,:), position(3*i-2,:))
   title('Theta 1 - Time vs Displacement','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta1 Joint Displacement (deg)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Velocity Graph
   hold on
   figure(5)
   plot(time(3*i-2,:), velocity(3*i-2,:))
   title('Theta 1 - Time vs Velocity','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta1 Joint Velocity (deg/s)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Acceleration Graph
   hold on
   figure(6)
   plot(time(3*i-2,:), acceleration(3*i-2,:))
   title('Theta 1 - Time vs Acceleration','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta1 Acceleration (deg/s^2)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   %%   Theta 2 Plots
   % Plotting Displacement Graph
   hold on
   figure(7)
   plot(time(3*i-1,:), position(3*i-1,:))
   title('Theta 2 - Time vs Displacement','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta2 Joint Displacement (deg)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Velocity Graph
   hold on
   figure(8)
   plot(time(3*i-1,:), velocity(3*i-1,:))
   title('Theta 2 - Time vs Velocity','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta2 Joint Velocity (deg/s)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
   
   % Plotting Acceleration Graph
   hold on
   figure(9)
   plot(time(3*i-1,:), acceleration(3*i-1,:))
   title('Theta 2 - Time vs Acceleration','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Theta2 Acceleration (deg/s^2)','fontsize', 12)
   legend('Joint 1', 'Joint 2', 'Joint 3', 'Joint 4', 'Joint 5', 'Joint 6')
   axis square
end

%%  Comet3D Plotting of e.e using link 1
%   Theta 1 is position(1,:)
%   Theta 2 is position(2,:)
%   Prismatic Length is position(3,:)


b0bi = L*[cosd(Beta(1)),cosd(Beta(2)),cosd(Beta(3)),cosd(Beta(4)),cosd(Beta(5)),cosd(Beta(6));...
          sind(Beta(1)),sind(Beta(2)),sind(Beta(3)),sind(Beta(4)),sind(Beta(5)),sind(Beta(6));
          0,0,0,0,0,0];

sum = zeros(3, length(position));
for i = 1:(n/3)
    th1 = position(3*i-2,:);
    th2 = position(3*i-1,:);
    lengths = position(3*i,:);

    prismatic = [L.*sind(Betad(i)) + lengths.*(sind(Betad(i)).*cosd(th2)+cosd(Betad(i)).*cosd(th1).*sind(th2))
                 -L.*cosd(Betad(i)) - lengths.*(cosd(Betad(i)).*cosd(th2)-sind(Betad(i)).*cosd(th1).*sind(th2))
                 lengths.*sind(th1).*sind(th2)];
    pi = b0bi(:,i) + prismatic;
    sum = sum + pi;
end
sum = sum./6;

figure(10)
comet3(sum(1,:),sum(2,:),sum(3,:))


