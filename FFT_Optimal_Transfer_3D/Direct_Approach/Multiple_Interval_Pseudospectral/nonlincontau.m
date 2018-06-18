%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% nonlincon.m
% nonlinear constraint function for multiple-interval pseudospectral methods
% defect equality constraints
% continuity equality constraints
% boundary equality constraints
% path inequality constraints
%--------------------------------------------------------------------------
% [c,ceq] = nonlincon(x,p)
%   x: scaled optimization variable vector
%   p: parameter structure
%   c: inequality constraint vector
% ceq: equality constraint vector
%--------------------------------------------------------------------------

function [c,ceq] = nonlincontau(x,p)

global Umax Rmin;
    % initialize constraints
    c = []; ceq = []; contin = []; in_contin = [];
    % extract from optimization vector
    % The optimization vector is a column vector containing the starting
    % column as states, i.e., ns*nt (number of n odes per interval are nt)
    
    X = reshape(x(1:p.ns*p.nt),p.nt,[]); % states % are reshaped into a matrix containing the state nodes as the rows and the node intervals as the columns 
    U = reshape(x(p.ns*p.nt+1:end),p.nt,[]); % controls % Control starts after the state variables are defined and are arranged in the same way as above
    
    %--- equality constraints
    % determine defect constraints for each segment
    for i = 1:length(p.Narray)
        Xi = X(p.cumN(i)+1:p.cumN(i+1),:); % states for segment i
        Ui = U(p.cumN(i)+1:p.cumN(i+1),:); % controls for segment i
        D = p.D{i}; % differentiation matrix
        F = Fmatrix(p.tau{i},Xi,Ui,p.t0(i),p.tf(i),p); % differential matrix
        defect{i,1} = D*Xi - F; % defect constraints
    end
    ceq = [ceq;reshape(cell2mat(defect),[],1)];
    
    % state and control continuity constraints between segments
    for i = 1:length(p.Narray)-1
        contin{i,1} = [X(p.cumN(i+1),:) - X(p.cumN(i+1)+1,:),...
            U(p.cumN(i+1),:) - U(p.cumN(i+1)+1,:)];
    end
    ceq = [ceq;reshape(cell2mat(contin),[],1)];

    % boundary constraints
    boundary = eval([p.boundary,'(X,U,p)']);
    ceq = [ceq;reshape(boundary,[],1)];
    
    Rx = X(:,1);
    Ry = X(:,3);
    Rz = X(:,5);
    
    Ux = U(:,1);
    Uy = U(:,2);
    Uz = U(:,3);
    
    
    Rnorm = (Rx.^2 + Ry.^2 + Rz.^2).^(0.5);
    Unorm = (Ux.^2 + Uy.^2 + Uz.^2).^(0.5);
    c = [c;Unorm-Umax];
    c = [c;-Rnorm+Rmin];
    
end