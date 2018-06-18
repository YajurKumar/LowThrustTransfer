%% The Coplanar Case Boundary Condition Setting %%
% Change this file if you want to change the initial and final conditions for the operation
% x(1)=x position, x(2)=x velocity, x(3)=y position, x(4)=y velocity
% 'a' stands for initial conditions, 'b' stands for final conditions 
% It is assumed that costate L1(0) = 0 or xa(5)=0
% L1=xb(5); L2=xb(6); L3=xb(7); L4=xb(8)

function res = BVP_bc(xa,xb,tf)

global x0 Vx0 y0 Vy0 xf Vxf yf Vyf N

res = [xa(1)-x0
    xa(2)-Vx0
    xa(3)-y0
    xa(4)-Vy0
    %xa(5)-0
    xb(1)-xf
    xb(2)-Vxf
    xb(3)-yf
    xb(4)-Vyf
    xb(5)*xb(2)+3*N*N*xb(6)*xb(1)+2*N*xb(6)*xb(4)+xb(7)*xb(4)-2*N*xb(8)*xb(2)-0.5*(xb(6)^2)-0.5*(xb(8)^2)
    ];

