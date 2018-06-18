%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% MCW_mayer.m
% Mayer term
%--------------------------------------------------------------------------
% m = BD_mayer(x,p)
% x: scaled optimization vector
% p: parameter structure
% m: Mayer term
%--------------------------------------------------------------------------

function m = MCW_mayer(x,p)
%     % extract from optimization vector
%     X = reshape(x(1:p.ns*p.nt),p.nt,[]); % states
%     U = reshape(x(p.ns*p.nt+1:end),p.nt,[]); % controls
%     % Lagrange term values for this segment
%     
%     xstart = 1;
%     xend   = xstart + (p.ns/6)*p.nt - 1;
%     vxstart = xend + 1;
%     vxend   = vxstart + (p.ns/6)*p.nt - 1;
%     ystart = vxend + 1;
%     yend   = ystart + (p.ns/6)*p.nt - 1;

%     vystart= yend + 1;
%     vyend  = vystart + (p.ns/6)*p.nt - 1;
%     zstart= vyend + 1;
%     zend  = zstart + (p.ns/6)*p.nt - 1;
%     vzstart= zend + 1;
%     vzend  = vzstart + (p.ns/6)*p.nt - 1;
%     
%     uxstart = p.ns*p.nt+1;
%     uxend   = uxstart + ((p.nu/3)*p.nt) - 1;
%     uystart = uxend + 1;
%     uyend   = uystart + ((p.nu/3)*p.nt) - 1;
%     uzstart = uyend + 1;
%     uzend   = uzstart + ((p.nu/3)*p.nt) - 1; 
%     
%     Rx = reshape(x(xstart:xend),p.nt,[]);
%     Vx = reshape(x(vxstart:vxend),p.nt,[]);
%     Ry = reshape(x(ystart:yend),p.nt,[]);
%     
%     Vy = reshape(x(vystart:vyend),p.nt,[]);
%     Rz = reshape(x(zstart:zend),p.nt,[]);
%     Vz = reshape(x(vzstart:vzend),p.nt,[]);
%     
%     Ux = reshape(x(uxstart:uxend),p.nt,[]);
%     Uy = reshape(x(uystart:uyend),p.nt,[]);
%     Uz = reshape(x(uzstart:uzend),p.nt,[]);

    %m = Rx(end)^2 + Ry(end)^2 + Rz(end)^2 + Vx(end)^2 + Vy(end)^2 + Vz(end)^2;
    m = 0; % none in this problem
end