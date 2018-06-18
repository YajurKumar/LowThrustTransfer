%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% objective.m
% objective function for multiple-interval pseudospectral methods
% calculates Lagrange and Mayer terms
%--------------------------------------------------------------------------

function fobj = objective(x,p)
    fobj = 0; % initialize
    % determine integral value
    for i = 1:length(p.Narray)
        h = eval([p.lagrange,'(x,p,i)']); % obtain segment Lagrange term values
        fobj = fobj + (p.tf(i)-p.t0(i))/2*p.w{i}'*h; % Gaussian quadrature  
    end
    % determine Mayer term value
    M = eval([p.mayer,'(x,p)']);
    fobj = fobj + M;
end