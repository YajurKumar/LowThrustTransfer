%% ODE For Hill-Clohessey-Wiltshire States and Co-States Equations: Coplanar Case, 
% Control Input applied along x and y axis
% States: x(1)=x, x(2)=Vx, x(3)=y, x(4)=Vy 
% Co-states: x(5)=L1, x(6)=L2, x(7)=L3, x(8)=L4

function dxdtau = BVP_ode(tau,x,tf)
format long;
global N
dxdt = [x(2)
    3*N*N*x(1)+2*N*x(4)-x(6)
    x(4)
    -2*N*x(2)-x(8)
    -3*N*N*x(6)
    -x(5)+2*N*x(8)
    0
    -2*N*x(6)-x(7)];
% multiplied with tf since, the time has been non-dimensionized
dxdtau = tf*dxdt;
return