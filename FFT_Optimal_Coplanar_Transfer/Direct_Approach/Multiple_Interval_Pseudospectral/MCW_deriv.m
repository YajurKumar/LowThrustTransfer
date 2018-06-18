%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% MCW_deriv.m
% Modified Clohessey-Wiltshire Derivative Function
% outputs a derivative function vector based on a single point in time
%--------------------------------------------------------------------------
% f = MCW_deriv(tau,X,U,t0,tf,p)
% tau: scaled nodes
%  X: state
%  U: control
% t0: initial time
% tf: final time
%  p: parameter structure
%  f: derivative value
%--------------------------------------------------------------------------

function f = MCW_deriv(tau,X,U,t0,tf,p)


global Re;

format long;

% J2 Value
J2    = 1.083e-3;

% Gravitational Parameter
mu    = 398600.4;

% Initial Radius of Reference Orbit
global r_0;

% Mean Motion
global n_0;

t = tf-t0;

% Inclination of the Orbit
i_d         = 51.6;
i           = deg2rad(i_d);

% Initial Argument of Latitude (Sum of Argument of Perigee and True
% Anomaly)
theta_0     = pi/6;
theta       = n_0*t + theta_0;


K = 6*J2*mu*(Re^2)/(r_0^5);

Xe = (3*J2*Re^2/(2*r_0))*((1/3)*((sin(i))^2)*((cos(theta))^2)+(1/3)*((sin(i))^2)-1+(1-(2/3)*(sin(i))^2)*cos(theta));

r     = r_0 + Xe;
r_dot = ((3/2)*J2*(Re^2)/(2*r_0))*((-1/3)*n_0*(sin(i)^2)*sin(2*theta)-n_0*(1-(2/3)*((sin(i))^2))*sin(theta));

% Computation of Angular Momentum
h_0     = sqrt(r_0*mu);
h_dot   = -(3/2)*J2*mu*Re^2*sin(2*theta)*((sin(i))^2)/(r^3);
h       = h_0 + ((3/4)*J2*mu*(Re^2)/(r_0*h_0))*(sin(i))^2*(cos(2*n_0*t)-1);

f_h       = -3/2*J2*mu*(Re^2)*sin(theta)*sin(2*i)/(r^4); 
gamma_dot = (r/h)*f_h;

omega_x = gamma_dot;
omega_y = 0;
omega_z = h/(r^2);

dummy = (n_0*h*(r^3)*cos(theta)-(h_dot*(r^3)+3*h*(r^2)*r_dot)*sin(theta))/((h^2)*(r^6));

omega_x_dot  = -3/2*J2*mu*Re^2*sin(2*i)*dummy;  
omega_y_dot  = 0;
omega_z_dot  = (h_dot*r-2*h*r_dot)/(r^3);

% Shorthands for State-Matrix
q_1 = omega_z^2 + 2*(n_0^2) + K*(1-3*(sin(i))^2*(sin(theta))^2);
q_2 = omega_z_dot + K*(sin(i))^2*sin(2*theta);
q_3 = -omega_x*omega_z + K*(sin(2*i)*sin(theta));
q_4 = -omega_z_dot + K*(sin(i))^2*sin(2*theta);
q_5 = omega_z^2 + omega_x^2 - n_0^2 + K*((-1/4)+(sin(i)^2)*((7/4)*(sin(theta))^2-(1/2)));
q_6 = omega_x_dot + K*(-1/4)*sin(2*i)*sin(theta);
q_7 = -omega_x*omega_z + K*sin(2*i)*sin(theta);
q_8 = -omega_x_dot + K*(-1/4)*sin(2*i)*cos(theta);
q_9 = omega_x^2 - n_0^2 + K*((-3/4)+(sin(i))^2*((5/4)*(sin(theta))^2+(1/2)));


% Initializing the differential equations
f(1,1) = X(2);
f(1,2) = q_1*X(1)+q_2*X(3)+2*omega_z*X(4)+q_3*X(5)+U(1);
f(1,3) = X(4);
f(1,4) = q_4*X(1)-2*omega_z*X(2)+q_5*X(3)+q_6*X(5)+2*omega_x*X(6)+U(2);
f(1,5) = X(6);
f(1,6) = q_7*X(1)+q_8*X(3)-2*omega_x*X(4)+q_9*X(5)+U(3);

end