%% Function Description
% This function is currently configured to use only symbols
% Outputs: one row of Jx, Jq 
% Ppx, Ppy, Pbx, Pby correspond to one joint at a time
% Ppx, Ppx = x & y distance from P0 -> Pi
% Pbx, Pby = x & y distance from B0 -> Bi 
% L = length of prismatic joint (Bi->Pi)
% Px, Py, Pz, alpha, beta, gamma = end effector positions & orientations (B0->P0)
% Note: Rotation matrix is currently configured to recieve angles in radians

function [Jx,Jq]=Jacobians_symbolic
%% Jacobian Jq
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
  
%% Jacobian Jx 
syms Px1 Py1 Pz1 Ppx1 Ppy1 Pbx1 Pby1 alpha1 beta1 gamma1 

% R0_ee=[cosd(alpha1)*cosd(beta1), cosd(alpha1)*sind(beta1)*sind(gamma1)-sind(alpha1)*cosd(gamma1), cosd(alpha1)*sind(beta1)*cosd(gamma1)+sind(alpha1)*sind(gamma1)
%            sind(alpha1)*cosd(beta1), sind(alpha1)*sind(beta1)*sind(gamma1)+cosd(alpha1)*cosd(gamma1), sind(alpha1)*sind(beta1)*cosd(gamma1)-cosd(alpha1)*sind(gamma1)
%            -sind(beta1)           , cosd(beta1)*sind(gamma1)                                    , cosd(beta1)*cosd(gamma1)];
R0_ee=[cos(alpha1)*cos(beta1), cos(alpha1)*sin(beta1)*sin(gamma1)-sin(alpha1)*cos(gamma1), cos(alpha1)*sin(beta1)*cos(gamma1)+sin(alpha1)*sin(gamma1)
           sin(alpha1)*cos(beta1), sin(alpha1)*sin(beta1)*sin(gamma1)+cos(alpha1)*cos(gamma1), sin(alpha1)*sin(beta1)*cos(gamma1)-cos(alpha1)*sin(gamma1)
           -sin(beta1)           , cos(beta1)*sin(gamma1)                                    , cos(beta1)*cos(gamma1)];
       
R0_p0= R0_ee; 
P_b0_p0= [Px1;Py1;Pz1];
P_p0_pi= [Ppx1;Ppy1;0]; 
P_b0_bi= [Pbx1;Pby1;0];
P_bi_pi=simplify(sum(expand((P_b0_p0+R0_p0*P_p0_pi-P_b0_bi).^2)));

%Gives you one row of the Jx matrix (ex. corresponding to joint 1)
Jx = [diff(P_bi_pi,Px1),diff(P_bi_pi,Py1),diff(P_bi_pi,Pz1),diff(P_bi_pi,alpha1),diff(P_bi_pi,beta1),diff(P_bi_pi,gamma1)];

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
return
