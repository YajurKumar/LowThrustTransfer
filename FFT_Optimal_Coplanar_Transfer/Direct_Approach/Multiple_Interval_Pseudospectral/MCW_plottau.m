%--------------------------------------------------------------------------
% MCW_plot.m
% MCW plotting function
%--------------------------------------------------------------------------
% MCW_plot(t,X,U,f,p)
% t: time
% X: state
% U: control
% f: objective function value
% p: parameter structure
%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 15, 2017
%
%--------------------------------------------------------------------------
function MCW_plottau(t,X,U,f,p)

global Umax ;
% interpolate the solution with the specified polynomials
interpN = 20000; % number of linearly spaced interpolation points
for i = 1:length(p.Narray)
    % interpolate based on method
    tarray1{i} = linspace(p.t0(i),p.tf(i),interpN);

    interpX11{1,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),1)',tarray1{i});
    interpX12{2,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),2)',tarray1{i});
    interpX13{3,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),3)',tarray1{i});
    interpX14{4,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),4)',tarray1{i});
    interpX15{5,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),5)',tarray1{i});
    interpX16{6,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),6)',tarray1{i});
    interpX17{6,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),7)',tarray1{i});
    interpX18{6,i} = LagrangeInter(p.t{i}',X(p.cumN(i)+1:p.cumN(i+1),8)',tarray1{i});
    
    interpU11{1,i} = LagrangeInter(p.t{i}',U(p.cumN(i)+1:p.cumN(i+1),1)',tarray1{i});
    interpU12{2,i} = LagrangeInter(p.t{i}',U(p.cumN(i)+1:p.cumN(i+1),2)',tarray1{i});
    interpU13{3,i} = LagrangeInter(p.t{i}',U(p.cumN(i)+1:p.cumN(i+1),3)',tarray1{i});

end

% create column vectors
tarray = cell2mat(tarray1)';

% interpX = cell2mat(interpX1)';
interpX1 = cell2mat(interpX11)';
interpX2 = cell2mat(interpX12)';
interpX3 = cell2mat(interpX13)';
interpX4 = cell2mat(interpX14)';
interpX5 = cell2mat(interpX15)';
interpX6 = cell2mat(interpX16)';
tmatrix  = cell2mat(interpX17)';
tfmatrix = cell2mat(interpX18)';

interpU1 = cell2mat(interpU11)';
interpU2 = cell2mat(interpU12)';
interpU3 = cell2mat(interpU13)';


figure
plot(tarray,tmatrix,'k')
ylabel('Time Taken (s)')
xlabel('Fraction of Specified Time (s)')
legend('x','y','z')
grid on;
title('Time taken trajectory')
set(gca,'fontsize',18)


figure
plot(tarray,tfmatrix,'k')
ylabel('Final Time Taken (s)')
xlabel('Fraction of Specified Time (s)')
legend('x','y','z')
grid on;
title('Final Time taken trajectory')
set(gca,'fontsize',18)


figure
plot(tarray,interpX1*Umax,'k-',tarray,interpX3*Umax,'k:',tarray,interpX5*Umax,'k-.','LineWidth',2)
ylabel('Relative Distance Components (km)')
xlabel('Time (s)')
legend('x','y','z')
grid on;
title('Relative Distance V. Time')
set(gca,'fontsize',18)

figure
plot(tarray,interpX2*Umax,'k-',tarray,interpX4*Umax,'k:',tarray,interpX6*Umax,'k-.','LineWidth',2)
grid on;
ylabel('Relative Velocity Components (km/s)')
xlabel('Time (s)')
legend('v_x','v_y','v_z')
title('Relative Velocity V. Time')
set(gca,'fontsize',18)

interpU = (interpU1.^2 + interpU2.^2 + interpU3.^2).^(0.5) ;
figure
plot(tarray,interpU1*Umax,'k-.',tarray,interpU2*Umax,'k:',tarray,interpU3*Umax,'k--',tarray,interpU*Umax,'k','LineWidth',2)
grid on;
xlabel('Time (s)');
ylabel('Control Input (km/s^2)');
title('Control Input V. Time')
legend('u_x','u_y','u_z','U');
set(gca,'fontsize',18)

figure
plot3(interpX1*Umax,interpX3*Umax,interpX5*Umax,'k','LineWidth',2)
grid on;
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Geometry of the Relative Orbit')
set(gca,'fontsize',18)

rho = (interpX1.^2 + interpX3.^2 + interpX5.^2).^(1/2);

figure
plot(tarray,rho*Umax,'k','LineWidth',2)
grid on;
xlabel('Time (s)')
ylabel('Relative Distance (km)')
title('Relative Distance Trajectory')
set(gca,'fontsize',18)


 end