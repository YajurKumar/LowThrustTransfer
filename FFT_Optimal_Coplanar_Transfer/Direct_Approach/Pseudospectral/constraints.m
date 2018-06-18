%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : May 5, 2017
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
%  g: derivative value
%--------------------------------------------------------------------------

function c = constraints(tau,X,U,t0,tf,p)

format long;
global Umax;
% Initializing the differential equations
c = (U(1).^2+U(2).^2+U(3).^2).^(0.5) - Umax;
end