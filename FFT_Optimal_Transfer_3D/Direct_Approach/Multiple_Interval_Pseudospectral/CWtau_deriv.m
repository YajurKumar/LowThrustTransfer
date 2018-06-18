% X(1)=x; X(2)=x'; X(3)=y; X(4)=y'; X(5)=z; X(6)=z'; X(7)=t; X(8)=tf;


function f = CWtau_deriv(tau,X,U,tau0,tauf,p)
format long;
global Umax;
global n_0;
N=n_0;
f(1,1) =  X(2)*X(8);
f(1,2) =  (3*N*N*X(1)+2*N*X(4))*X(8)+U(1);
f(1,3) =  X(4)*X(8);
f(1,4) =  -2*N*X(2)*X(8)+U(2);
f(1,5) =  X(6)*X(8);
f(1,6) =  -N*N*X(5)*X(8)+U(3);
f(1,7) =  X(8); % dt/dtau or dX(7)/dt
f(1,8) =  0; % dtf/dtau or dX(8)/dt
end
