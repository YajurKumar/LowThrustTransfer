%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% Fmatrix.m
% determines differential matrix for multiple-interval psuedospectral
% method
%--------------------------------------------------------------------------
% F = Fmatrix(tau,X,U,t0,tf,p)
% tau: nodes
%   X: states
%   U: control
%  t0: initial time
%  tf: final time
%   p: parameter structure
%   F: differential matrix
%--------------------------------------------------------------------------

function G = Gmatrix(tau,X,U,t0,tf,p)
    % number of nodes
    N = length(tau)-1;
    % initialize differential matrix
    f = zeros(N+1,size(X,2));
    % determine differential matrix
    for i = 1:N+1
        f(i,:) = eval([p.U_max,'(tau(i),X(i,:),U(i,:),t0,tf,p)']);
    end
    G = ((tf-t0)/2)*f; % scale
end