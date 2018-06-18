%% (C) Yajur Kumar, 2018. All rights reserved.
%v0.5. June 2018.
%% Program for optimizing the control (specific thrust) for relative coplanar transfer in the orbit
% The problem is to rendezvous the maneuvering spacecraft with the target 
% spacecraft while optimizing the overall fuel consumption. Since, the 
% transfer is coplanar, no out of plane terms are considered.  
% This is a free final time optimal control problem with boundary conditions specified. 

clc;
clear all;
close all;
format long;

% Declaring radius of target orbit (a), standard gravitational parameter
% for Earth (mu) and orbital frequency for target orbit (N) as global variables. 
global a mu N
a = 6858.14;
mu = 398600.4;
N = sqrt(mu/(a^3));

%----------------------------------------------------------------------------
%% Boundary Conditions
% x is in the radial direction directed outward from the center of Earth to the target
% y is along the velocity vector of the target
% Subscript 0 and f denotes initial condition at time t=0 and final
% condition at time t=t_f, respectively.
% Prefix V specifies the velocity
%----------------------------------------------------------------------------

global x0 Vx0 y0 Vy0 xf Vxf yf Vyf

% Initial Conditions
x0  = 10;   % in km
Vx0 = 0.01; % in km/s
y0  = 0.01;   % in km
Vy0 = 0.001;    % in km/s

% Final Conditions
xf  = 0;    % in km
Vxf = 0;    % in km/s
yf  = 0;    % in km
Vyf = 0;    % in km/s

% initial time
t0 = 0;

%----------------------------------------------------------------------------
%% Initial Guesses
%----------------------------------------------------------------------------

% initial guess for state and costate variables
init = [x0 Vx0 y0 Vy0 0 0 0 0];

% initial guess for mission time (s)
tf_guess = 500; 

Nt = 10; % number of data points

tau = linspace(0,1,Nt)'; % nondimensionalized time vector

solinit = bvpinit(tau,init,tf_guess);


options = bvpset('Stats','on','RelTol',1e-6,'NMax',1e6);

tic;
sol = bvp4c(@BVP_ode,@BVP_bc,solinit,options);
toc;

tf = sol.parameters(1);
Z = deval(sol,tau);


% Convert back to dimensional time for plotting
t = t0 - tau.*(tf-t0);
mission_time = t(end)

% Extract the solution for each state variable from the matrix Z:
x1 = Z(1,:);
x2 = Z(2,:);
x3 = Z(3,:);
x4 = Z(4,:);
u_x = -Z(6,:);
u_y = -Z(8,:);



figure
subplot(3,2,1);
plot(t,u_x)
set(gca,'fontsize',14)
ylabel('u_x (km/s^2)')
xlabel('t (s)')

subplot(3,2,2);
plot(t,u_y)
set(gca,'fontsize',14)
ylabel('u_y (km/s^2)')
xlabel('t (s)')

subplot(3,2,3);
plot(t,x1,t,x3)
set(gca,'fontsize',14)
legend('x','y')
ylabel('x,y (km)');
xlabel('t (s)');

r = (x1.^2 + x3.^2).^(1/2);
subplot(3,2,4);
plot(t,r)
set(gca,'fontsize',14)
ylabel('Relative distance (km)');
xlabel('t (s)');

subplot(3,2,5);
plot(t,x2,t,x4) % Trajectory of relative velocities
set(gca,'fontsize',14)
ylabel('v_x,v_y (km/s)');
xlabel('t (s)');

subplot(3,2,6);
plot(x1,x3)
set(gca,'fontsize',14) % Geometry of the maneuvering spacecraft orbit
xlabel('x (km)');
ylabel('y (km)');

