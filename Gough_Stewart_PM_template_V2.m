%1.1. Define size of figure and create figure handle (DO NOT MODIFY)
close all, clear all;

set(0,'Units','pixels');
dim = get(0,'ScreenSize'); 
fig_handle = figure('doublebuffer','on','Position',[0,35,dim(3),dim(4)-100],...
    'Name','3D Object','NumberTitle','off');
set(gcf,'color', [1 1 1]) %Background Colour

%1.2 Define the light in the figure (CHANGE POSITION VECTOR IF FIGURE IS TOO BRIGHT/DARK)
set(fig_handle,'Renderer','zbuffer','doublebuffer','off')
light('color',[.5,.5,.5],'position',[0,1,3],'Style','infinite')
lighting gouraud
daspect([1 1 1]);
axis off
view(60,30)

%Arrows (CHANGE PARAMETERS IF THEY ARE TOO SMALL OR TOO BIG)
% You need to have the file arrow3 in the same directory
   arrow_length=400; hold on
   line([0,0],[0,0], [0,arrow_length]); text(0,0,arrow_length*1.1,'z_0','FontSize',14); 
   line([0,0],[0,arrow_length],[0,0]); text(0,arrow_length*1.1, 0,'y_0','FontSize',14); 
   line([0,arrow_length],[0,0],[0,0]); text(arrow_length*1.1, 0, 0,'x_0','FontSize',14); 

%% Convert figure into Object (LOAD YOUR PARTS)       
load('Platform.mat');
setappdata(0,'object_data',object);
object = getappdata(0,'object_data');
obj{1}=object; 

load('Joints.mat');
setappdata(0,'object_data',object);
object = getappdata(0,'object_data');
obj{2}=object; 

load('Base.mat');
setappdata(0,'object_data',object);
object = getappdata(0,'object_data');
obj{3}=object; 

load('PrismaticJointRod.mat');
setappdata(0,'object_data',object);
object = getappdata(0,'object_data');
obj{4}=object;
obj{5}=object;
obj{6}=object;
obj{7}=object;
obj{8}=object;
obj{9}=object;

load('PrismaticJointShaft.mat');
setappdata(0,'object_data',object);
object = getappdata(0,'object_data');
obj{10}=object;
obj{11}=object;
obj{12}=object;
obj{13}=object;
obj{14}=object;
obj{15}=object;

for i=1:15  %(CHANGE if you have 10 parts, change it to i=1:10)
    q(i) = patch('faces', obj{i}.F, 'vertices', obj{i}.V);
    set(q(i),'EdgeColor','none');
end

%Set colour to the componenets (CHANGE colours of new parts)
set(q(1),'FaceColor', [1,0.242,0.293]);
set(q(2),'FaceColor', [.6,0.6,1]);
set(q(3),'FaceColor', [0.2,0.7,0.6]);
set(q(4:9),'FaceColor', [0.4,0.9,1]);
set(q(10:15),'FaceColor', [1,0.9,0.4]);


%%  KINEMATICS
%Characteristics of the manipulator
alpha=[-90,30,150]*pi/180;  % Definition of base platform
L=1500;                     % [mm] base length (centre to joint)
LL=1000;                      % mobile length (centre to joint)
P=3;                        % Link first length
PP =2.5;                    % Length second Length

%%%%%%%%%%%%%%%%%%%%
%INVERSE KINEMATICS%
%%%%%%%%%%%%%%%%%%%%
x0=500; y0=-1000; z0=2000;  %Random Position and Orientation of Mobile Platform
gamma=20; beta=30; alpha=-10;
R0_ee=[    cosd(beta)*cosd(gamma)                                 ,       -sind(beta)      , cosd(beta)*sind(gamma);...
       cosd(alpha)*sind(beta)*cosd(gamma)+sind(alpha)*sind(gamma) , cosd(alpha)*cosd(beta) , cosd(alpha)*sind(beta)*sind(gamma)-sind(alpha)*cosd(gamma);...
       sind(alpha)*sind(beta)*cosd(gamma)-cosd(alpha)*sind(gamma) , sind(alpha)*cosd(beta) , sind(alpha)*sind(beta)*sind(gamma)+cosd(alpha)*cosd(gamma)];

%WRITE HERE THE INVERSE KINEMATICS FOR ALL THE THREE BRANCHES, SOLVE FOR
%THE THREE ANGLES

%% My Stuff for IK
% Base Angles
Beta = [225
         315
         345
         75
         105
         195];
     
% Beta = angle to get to X axis for zero displacement config      
Betad = Beta+ 90;

%Platform Angles
anglePlat = [255
             285
             15
             45
             135
             165];

P = zeros(3,6);
theta1 = zeros(2,6);
Lnew = zeros(2,6);
theta2 = zeros(2,6);
R = zeros(3,3);
di = zeros(3,6);
unit = zeros(3,6);
Lnorm = zeros(1,6);         
b0bi = L*[cosd(Beta(1)),cosd(Beta(2)),cosd(Beta(3)),cosd(Beta(4)),cosd(Beta(5)),cosd(Beta(6));...
          sind(Beta(1)),sind(Beta(2)),sind(Beta(3)),sind(Beta(4)),sind(Beta(5)),sind(Beta(6));
          0,0,0,0,0,0];
p0pi = LL*[cosd(anglePlat(1)),cosd(anglePlat(2)),cosd(anglePlat(3)),cosd(anglePlat(4)),cosd(anglePlat(5)),cosd(anglePlat(6));...
          sind(anglePlat(1)),sind(anglePlat(2)),sind(anglePlat(3)),sind(anglePlat(4)),sind(anglePlat(5)),sind(anglePlat(6));
          0,0,0,0,0,0];

    


for i = 1:6
    di(:,i) = [x0; y0; z0] + LL*R0_ee*[cosd(anglePlat(i));sind(anglePlat(i));0] - b0bi(:,i); 
    
    Lnorm(i) = vecnorm(di(:,i));
    unit(:,i) = di(:,i)./Lnorm(i);

    %   Generate the position vector of spherical joint for IK
    P(:,i) = [x0; y0; z0] + LL*R0_ee*[cosd(anglePlat(i));sind(anglePlat(i));0];
    
    %   Theta 1 Inverse Kinematics equations for both solutions (solution 1
    %   in row 1, solution 2 in row 2)
    theta1(2,i) = atan2d(-P(3,i),-P(1,i)*cosd(Betad(i))-P(2,i)*sind(Betad(i)));
    theta1(1,i) = atan2d(P(3,i),P(1,i)*cosd(Betad(i))+P(2,i)*sind(Betad(i)));
    
    %   New prismatic joint length from IK 
    Lnew(1,i) = sqrt((P(1,i)*sind(Betad(i))-P(2,i)*cosd(Betad(i))-L)^2 + (P(3,i)*sind(theta1(1,i)) + P(1,i)*cosd(Betad(i))*cosd(theta1(1,i)) + P(2,i)*sind(Betad(i))*cosd(theta1(1,i)))^2);
    Lnew(2,i) = sqrt((P(1,i)*sind(Betad(i))-P(2,i)*cosd(Betad(i))-L)^2 + (P(3,i)*sind(theta1(2,i)) + P(1,i)*cosd(Betad(i))*cosd(theta1(2,i)) + P(2,i)*sind(Betad(i))*cosd(theta1(2,i)))^2);
    
    %   Theta 2 Inverse Kinematic Equations used as input for the FK
    theta2(1,i) = atan2d((P(3,i)*sind(theta1(1,i))+P(1,i)*cosd(Betad(i))*cosd(theta1(1,i))+P(2,i)*sind(Betad(i))*cosd(theta1(1,i)))/Lnew(1,i),(P(1,i)*sind(Betad(i))-P(2,i)*cosd(Betad(i))-L)/Lnew(1,i));
    theta2(2,i) = atan2d((P(3,i)*sind(theta1(2,i))+P(1,i)*cosd(Betad(i))*cosd(theta1(2,i))+P(2,i)*sind(Betad(i))*cosd(theta1(2,i)))/Lnew(2,i),(P(1,i)*sind(Betad(i))-P(2,i)*cosd(Betad(i))-L)/Lnew(2,i));
       

%Rotation Matrix to go from P to B0 
    R(:,:,i) = [cosd(Betad(i))*cosd(theta1(1,i))*cosd(theta2(1,i))-sind(Betad(i))*sind(theta2(1,i)),   -cosd(Betad(i))*sind(theta1(1,i)), sind(Betad(i))*cosd(theta2(1,i))+cosd(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i))
            cosd(Betad(i))*sind(theta2(1,i)) + sind(Betad(i))*cosd(theta1(1,i))*cosd(theta2(1,i)), -sind(Betad(i))*sind(theta1(1,i)), -cosd(Betad(i))*cosd(theta2(1,i))+sind(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i))
                                                              cosd(theta2(1,i))*sind(theta1(1,i)),  cosd(theta1(1,i)),                sind(theta1(1,i))*sind(theta2(1,i))];                                                                                       

% Rotation Matrix to go from P1 to B0                                                                                                                                   
%  R(:,:,i) = [cosd(Betad(i))*cosd(theta1(1,i))*cosd(theta2(1,i))-sind(Betad(i))*sind(theta2(1,i)), -sind(Betad(i))*cosd(theta2(1,i))-cosd(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i)), -cosd(Betad(i))*sind(theta1(1,i))
%              cosd(Betad(i))*sind(theta2(1,i)) + sind(Betad(i))*cosd(theta1(1,i))*cosd(theta2(1,i)), cosd(Betad(i))*cosd(theta2(1,i))-sind(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i)), -sind(Betad(i))*sind(theta1(1,i))
%                                                                cosd(theta2(1,i))*sind(theta1(1,i)),                                                sind(theta1(1,i))*sind(theta2(1,i)),cosd(theta1(1,i))];                                                                                       

end

%%%%%%%%%%%%%%%%%%%%
%FORWARD KINEMATICS%
FK = zeros(3,6);
for i=1:6    
    FK(1,i) =  L*sind(Betad(i))+ Lnorm(i)*(sind(Betad(i))*cosd(theta2(1,i)) + cosd(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i)));
    FK(2,i) = -L*cosd(Betad(i))- Lnorm(i)*(cosd(Betad(i))*cosd(theta2(1,i)) - sind(Betad(i))*cosd(theta1(1,i))*sind(theta2(1,i)));
    FK(3,i) = Lnorm(i)*sind(theta1(1,i))*sind(theta2(1,i));
end
%%%%%%%%%%%%%%%%%%%%

%%  Moving Parts
%Rename of all the vertices. This step is redundant obj{}.V will not longer be used.  
for i=1:15 %(CHANGE n=2 for the number of parts that you have)
    V{i} = obj{i}.V'; 
end

%Position of End Effector
newV{1} = R0_ee*V{1};
newV{1} = newV{1} + repmat([x0,y0,z0]',[1 length(V{1}(1,:))]); %Find new position of moving platform
newV{2} = R0_ee*V{2};
newV{2} = newV{2} + repmat([x0,y0,z0]',[1 length(V{2}(1,:))]); %Find new position of moving platform     


% NOTE IMPORTANT: For parts that rotate you have to multiply them by a
% rotation matrix first and then translate. For example, let T01a be the 
% rotation matrix that defines the frame 1 of a branch a with respect to the
% global reference frame located at the centre of the fixed platform. 
% This matrix involves the rotation of the first joint  

% Base Platform
newV{3} = V{3};

% Rotation of the Prismatic Joint
newV{4} = R(:,:,1)*V{4};
newV{5} = R(:,:,2)*V{5};
newV{6} = R(:,:,3)*V{6};
newV{7} = R(:,:,4)*V{7};
newV{8} = R(:,:,5)*V{8};
newV{9} = R(:,:,6)*V{9};
newV{10} = R(:,:,1)*V{10};
newV{11} = R(:,:,2)*V{11};
newV{12} = R(:,:,3)*V{12};
newV{13} = R(:,:,4)*V{13};
newV{14} = R(:,:,5)*V{14};
newV{15} = R(:,:,6)*V{15};

% Prismatic Joint Rods Translation
newV{4} = newV{4} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(255); LL*sind(255); 0],[1 length(V{4}(1,:))]);
newV{5} = newV{5} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(285); LL*sind(285); 0],[1 length(V{5}(1,:))]);
newV{6} = newV{6} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(15); LL*sind(15); 0],[1 length(V{6}(1,:))]);
newV{7} = newV{7} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(45); LL*sind(45); 0],[1 length(V{7}(1,:))]);
newV{8} = newV{8} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(135); LL*sind(135); 0],[1 length(V{8}(1,:))]);
newV{9} = newV{9} + repmat([x0; y0; z0] + R0_ee*[LL*cosd(165); LL*sind(165); 0],[1 length(V{9}(1,:))]);

% Prismatic Joint Shafts Translation
newV{10} = newV{10} + repmat([L*cosd(225), L*sind(225), 0]',[1 length(V{10}(1,:))]);
newV{11} = newV{11} + repmat([L*cosd(315), L*sind(315), 0]',[1 length(V{11}(1,:))]);
newV{12} = newV{12} + repmat([L*cosd(345), L*sind(345), 0]',[1 length(V{12}(1,:))]);
newV{13} = newV{13} + repmat([L*cosd(75),  L*sind(75), 0]',[1 length(V{13}(1,:))]);
newV{14} = newV{14} + repmat([L*cosd(105), L*sind(105), 0]',[1 length(V{14}(1,:))]);
newV{15} = newV{15} + repmat([L*cosd(195), L*sind(195), 0]',[1 length(V{15}(1,:))]); 

%newV{3} = newV{3} + repmat([x0,y0,z0]',[1 length(V{3}(1,:))]);
%newV{3} = T01a(1:3,1:3)*V{3}; %Find new orientation of link 1 branch 1
%newV{3} = newV{3} + repmat(T01a(1:3,4),[1 length(newV{3}(1,:))]); 


for ii=[1:15] %(CHANGE n=2 to the number of parts that you have)
    set(q(ii),'Vertices',newV{ii}(1:3,:)'); %Set the new position in the handle (graphical link)
end
             
drawnow
