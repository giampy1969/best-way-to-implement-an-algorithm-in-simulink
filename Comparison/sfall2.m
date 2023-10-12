function sfall2(block)

% Level-2 MATLAB file S-Function for a 2-state EKF for attitude estimation
% For more info on the algorithm, its implementation and usage, go to 
% the root folder and see the files README.txt, Presentation.pdf, and 
% ekf2sim_example.slx (read the help of the ekf2sim block therein). 
% Copyright 2014 The MathWorks, Inc.

setup(block);

%endfunction

function setup(block)

block.NumDialogPrms  = 7;

%% Register number of input and output ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 1;

%% Setup functional port properties to dynamically
%% inherited.
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

block.InputPort(1).Dimensions        = 3;
block.InputPort(1).DirectFeedthrough = false;

block.InputPort(2).Dimensions        = 3;
block.InputPort(2).DirectFeedthrough = false;

block.InputPort(3).Dimensions        = 3;
block.InputPort(3).DirectFeedthrough = true;

block.OutputPort(1).Dimensions       = 3;

%% Set block sample time to inherited
block.SampleTimes = [block.DialogPrm(6).Data 0];

%% Register methods
block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions',    @InitConditions);
block.RegBlockMethod('Outputs',                 @Output);
block.RegBlockMethod('Update',                  @Update);

%endfunction

function DoPostPropSetup(block)

%% Setup Dwork
block.NumDworks = 3;

block.Dwork(1).Name = 'P';
block.Dwork(1).Dimensions      = 4;
block.Dwork(1).DatatypeID      = 0;
block.Dwork(1).Complexity      = 'Real';
block.Dwork(1).UsedAsDiscState = true;

block.Dwork(2).Name = 'x';
block.Dwork(2).Dimensions      = 2;
block.Dwork(2).DatatypeID      = 0;
block.Dwork(2).Complexity      = 'Real';
block.Dwork(2).UsedAsDiscState = true;

block.Dwork(3).Name = 'vold';
block.Dwork(3).Dimensions      = 3;
block.Dwork(3).DatatypeID      = 0;
block.Dwork(3).Complexity      = 'Real';
block.Dwork(3).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

% define dimensions
nx=2;nu=3;ny=2;

block.Dwork(1).Data=reshape(block.DialogPrm(3).Data,nx*nx,1);
block.Dwork(2).Data=block.DialogPrm(4).Data;
block.Dwork(3).Data=block.DialogPrm(7).Data;

%endfunction

function Output(block)


% get inputs
ve = block.InputPort(3).Data;

% calculate heading
psi=atan2(ve(2),ve(1));

%% attach psi to x and sent all to output
block.OutputPort(1).Data = [block.Dwork(2).Data;psi];


function Update(block)

%% function [Pp,xp,psi]=ekf2tx(P,x,vold,W,V,T,G, Plim,ve,u,ab)

% define dimensions
nx=2;nu=3;ny=2;

% get states
P = reshape(block.Dwork(1).Data,2,2);
x = block.Dwork(2).Data;
vold = block.Dwork(3).Data;

% get parameters
W = block.DialogPrm(1).Data;
V = block.DialogPrm(2).Data;
T = block.DialogPrm(6).Data;
Plim=block.DialogPrm(5).Data;

% get inputs
u = block.InputPort(1).Data;
ab = block.InputPort(2).Data;
ve = block.InputPort(3).Data;

% calculate gps acceleration
agps = (ve-vold)/T;

% calculate heading
psi=atan2(ve(2),ve(1));

% rotate gps acceleration
rx=-cos(psi)*agps(1)-sin(psi)*agps(2);
ry=sin(psi)*agps(1)-cos(psi)*agps(2);
rz=9.80665-agps(3);

% solve for theta and phi given imu and gps accelerations
sth=(rx*ab(1)+rz*(abs(rx^2+rz^2-ab(1)^2)^0.5))/(realmin+rx^2+rz^2);
theta=atan2(sth*rx-ab(1),sth*rz);
rth=rx*sin(theta)+rz*cos(theta);
sph=(ry*ab(2)+rth*(abs(ry^2+rth^2-ab(2)^2)^0.5))/(realmin+ry^2+rth^2);
phi=atan2(-sph*ry+ab(2),sph*rth);

% this works as the "measured" output in the ekf
y = [phi;theta];

% states (note, phi and theta are re-initialized with the previous estimates)
phi=x(1);
theta=x(2);

% inputs
p=u(1);
q=u(2);
r=u(3);

% omega to attitude dot
Taw=[1 sin(phi)*tan(theta) cos(phi)*tan(theta); 0 cos(phi) -sin(phi)];

% x propagation
xdot=Taw*[p;q;r];
xm=x+T*xdot;

% A matrix
a = [ cos(phi)*tan(theta)*q-sin(phi)*tan(theta)*r      sin(phi)*(1+tan(theta)^2)*q+cos(phi)*(1+tan(theta)^2)*r
      -sin(phi)*q-cos(phi)*r                           0                                                      ];

% euler forward step
A = eye(2)+T*a; 

% find K
Pm=W+A*P*A';
VPm=V+Pm;
VPm_det=VPm(1,1)*VPm(2,2)-VPm(1,2)*VPm(2,1);
iVPm=[VPm(2,2) -VPm(1,2); -VPm(2,1) VPm(1,1)]/VPm_det;
K=Pm*iVPm;

% EKF update equations
Pp=(eye(2)-K)*Pm;
xp=xm+K*(y-xm);

% ball norm limiter
nrm=sqrt(Pp(:)'*Pp(:));
Pp=Pp*min(nrm,Plim)/max(nrm,realmin);

% store
block.Dwork(1).Data = reshape(Pp,nx*nx,1);
block.Dwork(2).Data = xp;
block.Dwork(3).Data = ve;

%endfunction

