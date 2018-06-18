%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% CW_deriv.m
% Clohessey-Wiltshire Derivative Function
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

function f = CW_deriv(tau,X,U,t0,tf,p)
format long;
global n_0
N=n_0;
f(1,1) =  X(2);
f(1,2) =  3*N*N*X(1)+2*N*X(4)+U(1);
f(1,3) =  X(4);
f(1,4) =  -2*N*X(2)+U(2);
f(1,5) =  X(6);
f(1,6) =  -N*N*X(5)+U(3);
end
