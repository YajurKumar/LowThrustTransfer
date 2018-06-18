%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% MCW_boundary.m
% Modified Clohessey-Wiltshire boundary constraint function
%--------------------------------------------------------------------------
% b = MCW_boundary(X,U,p)
% X: State
% U: Control
% p: Parameter structure
% b: Boundary constraint vector
%--------------------------------------------------------------------------

function b = MCW_boundary(X,U,p)
    % initialize boundary constraint
    b = [];
    % X1 state boundary constraints
    b = [b;X(1,1)-p.X1_0;X(end,1)-p.X1_f];
    % X2 state boundary constraints
    b = [b;X(1,2)-p.X2_0;X(end,2)-p.X2_f];
    % X3 state boundary constraints
    b = [b;X(1,3)-p.X3_0;X(end,3)-p.X3_f];
    % X4 state boundary constraints
    b = [b;X(1,4)-p.X4_0;X(end,4)-p.X4_f];
    % X5 state boundary constraints
    b = [b;X(1,5)-p.X5_0;X(end,5)-p.X5_f];
    % X6 state boundary constraints
    b = [b;X(1,6)-p.X6_0;X(end,6)-p.X6_f];
end