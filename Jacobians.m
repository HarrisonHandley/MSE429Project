%% Function Description
% Find: Jx, Jq, qd (joint rates), T (torque)
% Inputs Ppx, Ppy, Pbx, Pby correspond to one joint at a time
% Ppx, Ppx = x & y distance from P0 -> Pi
% Pbx, Pby = x & y distance from B0 -> Bi 
% L = length of prismatic joint (Bi->Pi)
% Px, Py, Pz, alpha, beta, gamma = end effector positions & orientations (B0->P0). These
% are all contained in the position matrix and can be subsituted in from the position matrix 
% T is a function of xd (task space velocities) which can be found under
% Path Generation 
% NOTE: Rotation matrix is currently configured to recieve angles in RADIANS

%% PATH GENERATION 
%Using this section to find the EE velocity (xd) and positions  
% xd contains the velocities for Px, Py,Pz,alpha,beta,gamma for every step
% dt=0.1; %stepsize
% P_ee = [0 0.05 0 0; 
%         0 0.05 0 0;
%         0.2 0.2 0.2 0.2;
%         0 45 45 0; 
%         0 0 0 0;
%         0 0.7068766 0.6332458 0];
% tf = [0.40385,0.8077,0.8077];
% 
% [position,xd,~,~]=via_points_match_VA(P_ee, tf, dt, 'prescribed',[0,0]);

%% Constants
L=1500;                     % [mm] base length (centre to joint)
LL=1000; 
% Base Angles
Beta = [225
         315
         345
         75
         105
         195];

%Platform Angles
anglePlat = [255
             285
             15
             45
             135
             165];
b0bi = L*[cosd(Beta(1)),cosd(Beta(2)),cosd(Beta(3)),cosd(Beta(4)),cosd(Beta(5)),cosd(Beta(6));...
          sind(Beta(1)),sind(Beta(2)),sind(Beta(3)),sind(Beta(4)),sind(Beta(5)),sind(Beta(6));
          0,0,0,0,0,0];
p0pi = LL*[cosd(anglePlat(1)),cosd(anglePlat(2)),cosd(anglePlat(3)),cosd(anglePlat(4)),cosd(anglePlat(5)),cosd(anglePlat(6));...
          sind(anglePlat(1)),sind(anglePlat(2)),sind(anglePlat(3)),sind(anglePlat(4)),sind(anglePlat(5)),sind(anglePlat(6));
          0,0,0,0,0,0];
      
%% Jacobian Jq
% L_sqr = Lnorm^2; 
% Derivative of Lnorm^2 gives you 2*Lnorm  
% This matrix is 6x6 containing diagonals with lengths L1,L2,L3 etc       
syms L1 L2 L3 L4 L5 L6 
Jq = 2*[L1,0,0,0,0,0;
      0,L2,0,0,0,0;
      0,0,L3,0,0,0;
      0,0,0,L4,0,0;
      0,0,0,0,L5,0;
      0,0,0,0,0,L6];
  
%% Jacobian Jx
Jx = zeros(6,6);
Jx = sym(Jx);
for i=1:6
[Jx(i,:)]=findJx(p0pi(1,i), p0pi(2,i), b0bi(1,i), b0bi(2,i));
end

%% Forward Velocity Problem (joint rate velocities as function of end effector velocity)
%xd contains the values for end effector velocity 
%xd = derivative([x,y,z,a,b,g]')
syms Dx Dy Dz Dalpha Dbeta Dgamma 
xd = [Dx 
      Dy 
      Dz 
      Dalpha 
      Dbeta 
      Dgamma];

% symbolic joint rate velocity matrix 
%qd = derivative([L1,L2,L3,L4,L5,L6]')
% calculated values for joint rate velocities 
qd = Jq^(-1)*(Jx*xd); 
qd = simplify(qd);


%% Forward Static Force Problem (required torques based on external force at end effector)
%Weight of mobile platform & things on top = external force
weight =271.4e-3; 
F = [0;0;weight;0;0;0];
T = transpose(Jq^(-1)*Jx)*F;

%% Function to populate Jx 
% Takes in platform and base positions for specific joint
% Returns: Jacobian Row for that joint 

function [Jx]=findJx(Ppx, Ppy, Pbx, Pby)
%% Jacobian Jx
syms Px1 Py1 Pz1 Ppx1 Ppy1 Pbx1 Pby1 a b g
R0_ee=[cos(a)*cos(b), cos(a)*sin(b)*sin(g)-sin(a)*cos(g), cos(a)*sin(b)*cos(g)+sin(a)*sin(g)
           sin(a)*cos(b), sin(a)*sin(b)*sin(g)+cos(a)*cos(g), sin(a)*sin(b)*cos(g)-cos(a)*sin(g)
           -sin(b)           , cos(b)*sin(g)                                    , cos(b)*cos(g)];
       
R0_p0= R0_ee; 
P_b0_p0= [Px1;Py1;Pz1];
P_p0_pi= [Ppx1;Ppy1;0]; 
P_b0_bi= [Pbx1;Pby1;0];
P_bi_pi=simplify(sum(expand((P_b0_p0+R0_p0*P_p0_pi-P_b0_bi).^2)));

%Gives you one row of the Jx matrix (ex. corresponding to joint 1)
Jx = [diff(P_bi_pi,Px1),diff(P_bi_pi,Py1),diff(P_bi_pi,Pz1),diff(P_bi_pi,a),diff(P_bi_pi,b),diff(P_bi_pi,g)];

%clear Px1 Py1 Pz1 Ppx1 Ppy1 Pbx1 Pby1 alpha1 beta1 gamma1 

Jx1 = subs(Jx,Ppx1,Ppx);
Jx1 = subs(Jx1,Ppy1,Ppy);
Jx1 = subs(Jx1,Pbx1,Pbx);
Jx1 = subs(Jx1,Pby1,Pby);
Jx = simplify(Jx1);
end





