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


%% TRAJECTORY GENERATION. Trajectory generation of all the joints
dt = 0.1; %stepsize
[position,velocity,acceleration,time] = via_points_match_VA(P_ee, tf, dt, 'prescribed',[0,0]);

%% Plotting Graphs
[n, j] = size(acceleration);    % rows, columns

for i = 1:(n/2)
   %%   Cartesian Coordinates Graphs
   % Plotting Displacement Graph
   hold on
   figure(1)
   plot(time(i,:), position(i,:)./10000)
   title('End Effector Cartesian Coordinates - Time vs Displacement','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Displacement (m)','fontsize', 12)
   legend('x', 'y', 'z')
   axis square
   
   % Plotting Velocity Graph
   hold on
   figure(2)
   plot(time(i,:), velocity(i,:)./10000)
   title('End Effector Cartesian Coordinates - Time vs Velocity','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Velocity (m/s)','fontsize', 12)
   legend('x', 'y', 'z')
   axis square
   
   % Plotting Acceleration Graph
   hold on
   figure(3)
   plot(time(i,:), acceleration(i,:)./10000)
   title('End Effector Cartesian Coordinates - Time vs Acceleration','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Acceleration (m/s^2)','fontsize', 12)
   legend('x', 'y', 'z')
   axis square
   
   %%   Angular Cartesian 
   % Plotting Displacement Graph
   hold on
   figure(4)
   plot(time(i+3,:), position(i+3,:))
   title('End Effector Angular Coordinates - Time vs Displacement','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Angular Displacement (deg)','fontsize', 12)
   legend('alpha', 'beta', 'gamma')
   axis square
   
   % Plotting Velocity Graph
   hold on
   figure(5)
   plot(time(i+3,:), velocity(i+3,:))
   title('End Effector Angular Coordinates - Time vs Velocity','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Angular Velocity (deg/s)','fontsize', 12)
   legend('alpha', 'beta', 'gamma')
   axis square
   
   % Plotting Acceleration Graph
   hold on
   figure(6)
   plot(time(i+3,:), acceleration(i+3,:))
   title('End Effector Angular Coordinates - Time vs Acceleration','fontsize', 14)
   xlabel('Time (s)','fontsize', 12)
   ylabel('Angular Acceleration (deg/s^2)','fontsize', 12)
   legend('alpha', 'beta', 'gamma')
   axis square
   
end

%%  Comet3D Plotting of e.e using link 1
%   Theta 1 is position(1,:)
%   Theta 2 is position(2,:)
%   Prismatic Length is position(3,:)

