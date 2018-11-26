%% Function Description
%Two options:
%1) use inputs (See line 192 in Gough_Stewart_PM_template_V2.m)
%2) use symbolic stuff
% This function is currently configured to use only symbols
% Outputs: Jx, Jq, symbolic matrix qd (joint rates) 

% For 2), Inputs Ppx, Ppy, Pbx, Pby correspond to one joint at a time
% Ppx, Ppx = x & y distance from P0 -> Pi
% Pbx, Pby = x & y distance from B0 -> Bi 
% Lnorm = length of prismatic joint (Bi->Pi)
% Px, Py, Pz, alpha, beta, gamma = end effector positions & orientations (B0->P0). These
% are all contained in the position matrix and need to be iterated through 
% Note: Rotation matrix is currently configured to recieve angles in
% radians

%% NEED TO DO!!
% 1. Make Jq a full matrix (6x6), currently 1x6
% 2. Find qd_sym (joint velocities) -> Forward Velocity Problem
% 3. Find T_sym (torques) -> Inverse Velocity Problem 

function [Jx,Jq,qd_sym,T_sym]=Fk_practice
%Uncomment to send actual inputs with values 
%function [Jx,Jq,qd,T_sym]=Fk_practice(Lnorm,Px, Py, Pz, Ppx, Ppy, Pbx, Pby, alpha, beta, gamma)

%% PATH GENERATION 
%Using this section to find the EE velocity (xd) 
% xd contains the velocities for Px, Py,Pz,alpha,beta,gamma for every step

dt=0.1; %stepsize
P_ee = [0 0.05 0 0; 
        0 0.05 0 0;
        0.2 0.2 0.2 0.2;
        0 45 45 0; 
        0 0 0 0;
        0 0.7068766 0.6332458 0];
tf = [0.40385,0.8077,0.8077];

[position,xd,~,~]=via_points_match_VA(P_ee, tf, dt, 'prescribed',[0,0]);


%% Jacobian Jq
% Left hand side (Li^2)
 %L_sqr = Lnorm^2; 
% Derivative of Lnorm^2 gives you 2*Lnorm  
% This matrix is 6x6 containing diagonals with lengths L1,L2,L3 etc 
syms L1 L2 L3 L4 L5 L6 
Jq = 2*[L1,0,0,0,0,0;
      0,L2,0,0,0,0;
      0,0,L3,0,0,0;
      0,0,0,L4,0,0;
      0,0,0,0,L5,0;
      0,0,0,0,0,L6];

% Uncomment to use actual values
% Jq = 2*[Lnorm(1),0,0,0,0,0;
%       0,Lnorm(2),0,0,0,0;
%       0,0,Lnorm(3),0,0,0;
%       0,0,0,Lnorm(4),0,0;
%       0,0,0,0,Lnorm(5),0;
%       0,0,0,0,0,Lnorm(6)];
  

%% Jacobian Jx 
syms Px1 Py1 Pz1 Ppx1 Ppy1 Pbx1 Pby1 alpha1 beta1 gamma1 

% R0_ee=[cosd(alpha1)*cosd(beta1), cosd(alpha1)*sind(beta1)*sind(gamma1)-sind(alpha1)*cosd(gamma1), cosd(alpha1)*sind(beta1)*cosd(gamma1)+sind(alpha1)*sind(gamma1)
%            sind(alpha1)*cosd(beta1), sind(alpha1)*sind(beta1)*sind(gamma1)+cosd(alpha1)*cosd(gamma1), sind(alpha1)*sind(beta1)*cosd(gamma1)-cosd(alpha1)*sind(gamma1)
%            -sind(beta1)           , cosd(beta1)*sind(gamma1)                                    , cosd(beta1)*cosd(gamma1)];
R0_ee=[cos(alpha1)*cos(beta1), cos(alpha1)*sin(beta1)*sin(gamma1)-sin(alpha1)*cos(gamma1), cos(alpha1)*sin(beta1)*cos(gamma1)+sin(alpha1)*sin(gamma1)
           sin(alpha1)*cos(beta1), sin(alpha1)*sin(beta1)*sin(gamma1)+cos(alpha1)*cos(gamma1), sin(alpha1)*sin(beta1)*cos(gamma1)-cos(alpha1)*sin(gamma1)
           -sin(beta1)           , cos(beta1)*sin(gamma1)                                    , cos(beta1)*cos(gamma1)];
       
R0_p0= R0_ee; 
P_b0_p0= [Px1,Py1,Pz1]';
P_p0_pi= [Ppx1,Ppy1,0]'; 
P_b0_bi= [Pbx1,Pby1,0]';
P_bi_pi=simplify(sum(expand((P_b0_p0+R0_p0*P_p0_pi-P_b0_bi).^2)));

%Gives you one row of the Jx matrix (ex. corresponding to joint 1)
Jx = [diff(P_bi_pi,Px1);diff(P_bi_pi,Py1);diff(P_bi_pi,Pz1);diff(P_bi_pi,alpha1);diff(P_bi_pi,beta1);diff(P_bi_pi,gamma1)]';

%clear Px1 Py1 Pz1 Ppx1 Ppy1 Pbx1 Pby1 alpha1 beta1 gamma1 

% Uncomment and Use Jx1 to calculate actual values instead of symbolic 
% Jx1 = subs(Jx,Px1,Px);
% Jx1 = subs(Jx1,Py1,Py);
% Jx1 = subs(Jx1,Pz1,Pz);
% Jx1 = subs(Jx1,Ppx1,Ppx);
% Jx1 = subs(Jx1,Ppy1,Ppy);
% Jx1 = subs(Jx1,Pbx1,Pbx);
% Jx1 = subs(Jx1,Pby1,Pby);
% Jx1 = subs(Jx1,alpha1,alpha);
% Jx1 = subs(Jx1,beta1,beta);
% Jx1 = subs(Jx1,gamma1,gamma);
% Jx1 = simplify(Jx1);

%% Forward Velocity Problem (joint rate velocities as function of end effector velocity)
%xd contains the values for end effector velocity 
%xd = [x.,y.,z.,a.,b.,g.]'
syms Dx Dy Dz Dalpha Dbeta Dgamma 
xd = [Dx Dy Dz Dalpha Dbeta Dgamma];
for i=1:length(xd(1,:))

% symbolic joint rate velocity matrix 
%qd = [L1.,L2.,L3.,L4.,L5.,L6.]
% This won't be right because Jx is only one row right now (1x6) instead of
% 6x6 
qd_sym = Jq^(-1)*(Jx*xd');

% calculated values for joint rate velocities 
% qd = Jq^(-1)*(Jx1*xd(:,i)); 
% qd = double(simplify(qd));

end 


%% Forward Static Problem (required torques based on external force at end effector)
%Weight of mobile platform & things on top = external force
syms weight % Change this
F = [0,0,weight,0,0,0]';
% This won't work right now because Jx is currently just one row (1x6)
% matrix and needs to be 6x6 
%T_sym = transpose(Jq^(-1)*Jx)*F;
T_sym = 0;



return
