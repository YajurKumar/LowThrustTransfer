%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% MCW_lagrange.m
% MCW Lagrange term
%--------------------------------------------------------------------------
% L = MCW_lagrange(x,p,i)
% x: scaled optimization vector
% p: parameter structure
% i: interval number
% L: Lagrange term
%--------------------------------------------------------------------------

function L = MCW_lagrange(x,p,i)
    % extract from optimization vector
    X = reshape(x(1:p.ns*p.nt),p.nt,[]); % states
    % U = reshape(x(p.ns*p.nt+1:end),p.nt,[]); % controls
    
    uxstart = p.ns*p.nt+1;
    uxend   = uxstart + ((p.nu/3)*p.nt) - 1;
    uystart = uxend + 1;
    uyend   = uystart + ((p.nu/3)*p.nt) - 1;
    uzstart = uyend + 1;
    uzend   = uzstart + ((p.nu/3)*p.nt) - 1; 
    
    Ux = reshape(x(uxstart:uxend),p.nt,[]);
    Uy = reshape(x(uystart:uyend),p.nt,[]);
    Uz = reshape(x(uzstart:uzend),p.nt,[]);
    
    % Lagrange term values for this segment
    L = 1/2*(Ux(p.cumN(i)+1:p.cumN(i+1),1).^2 + Uy(p.cumN(i)+1:p.cumN(i+1),1).^2 + Uz(p.cumN(i)+1:p.cumN(i+1),1).^2) ;
end