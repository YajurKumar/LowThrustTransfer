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
function MCW_plot(t,X,U,f,p)

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

interpU1 = cell2mat(interpU11)';
interpU2 = cell2mat(interpU12)';
interpU3 = cell2mat(interpU13)';

% figure
% plot(tarray,interpX)
% xlabel('Time (s)');
% ylabel('States');

% figure
% plot(tarray,interpX1)
% xlabel('Time (s)');
% ylabel('x (km)');
% title('Relative Distance along x axis V. Time')
% figure
% plot(tarray,interpX2)
% xlabel('Time (s)');
% ylabel('v_x (km/s)');
% title('Relative Velocity along x axis V. Time')
% figure
% plot(tarray,interpX3)
% xlabel('Time (s)');
% ylabel('y (km)');
% title('Relative Distance along y axis V. Time')
% figure
% plot(tarray,interpX4)
% xlabel('Time (s)');
% ylabel('v_y (km/s)');
% title('Relative Velocity along y axis V. Time')
% figure
% plot(tarray,interpX5)
% xlabel('Time (s)');
% ylabel('z (km)');
% title('Relative Distance along z axis V. Time')
% figure
% plot(tarray,interpX6)
% xlabel('Time (s)');
% ylabel('v_z (km/s)');
% title('Relative Velocity along z axis V. Time')
figure
plot(tarray,interpX1,'k-',tarray,interpX3,'k:',tarray,interpX5,'k-.','LineWidth',2)
ylabel('Relative Distance Components (km)')
xlabel('Time (s)')
legend('x','y','z')
grid on;
title('Relative Distance V. Time')
set(gca,'fontsize',18)

figure
plot(tarray,interpX2,'k-',tarray,interpX4,'k:',tarray,interpX6,'k-.','LineWidth',2)
grid on;
ylabel('Relative Velocity Components (km/s)')
xlabel('Time (s)')
legend('v_x','v_y','v_z')
title('Relative Velocity V. Time')
set(gca,'fontsize',18)

interpU = (interpU1.^2 + interpU2.^2 + interpU3.^2).^(0.5) ;
figure
plot(tarray,interpU1,'k-.',tarray,interpU2,'k:',tarray,interpU3,'k--',tarray,interpU,'k','LineWidth',2)
grid on;
xlabel('Time (s)');
ylabel('Control Input (km/s^2)');
title('Control Input V. Time')
legend('u_x','u_y','u_z','U');
set(gca,'fontsize',18)

figure
plot3(interpX1,interpX3,interpX5,'k','LineWidth',2)
grid on;
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')
title('Geometry of the Relative Orbit')
set(gca,'fontsize',18)

rho = (interpX1.^2 + interpX3.^2 + interpX5.^2).^(1/2);

figure
plot(tarray,rho,'k','LineWidth',2)
grid on;
xlabel('Time (s)')
ylabel('Relative Distance (km)')
title('Relative Distance Trajectory')
set(gca,'fontsize',18)


% figure
% [hAinterp,hd1interp,hd2interp] = plotyy(tarray,interpX1,tarray,interpX2,tarray,interpX3,tarray,interpX4,tarray,interpX5,tarray,interpX6,tarray,interpU)
% title('Plot of States with Legendre Interpolation')
% 
% figure
% [hA,hd1,hd2] = plotyy(t,X,t,U)
% title('Plot of States without Legendre Interpolation')
% set(hAinterp(1),'YTick',[])
% set(hAinterp(2),'YTick',[])
% set(hA(1),'YTick',[])
% set(hA(2),'YTick',[])
% 
% % interp style
% set(hd1interp(1),'linewidth',2,'color',xc);
% set(hd1interp(2),'linewidth',2,'color',vc);
% set(hd2interp(1),'linewidth',2,'color',uc);
% 
% % node style
% set(hd1(1),'linestyle','none','Marker','o','linewidth',2,'color',xc);
% set(hd1(2),'linestyle','none','Marker','o','linewidth',2,'color',vc);
% set(hd2(1),'linestyle','none','Marker','o','linewidth',2,'color',uc);
% 
% % y labels
% ylabel(hAinterp(1),'$x$ [m], $v$ [m/s]','interpreter','latex','fontsize',font1) % left y-axis
% ylabel(hAinterp(2),'$u$ [m/s$^2$]','interpreter','latex','fontsize',font1,'color',uc) % right y-axis
% 
% % 
% set(hAinterp(1),'ylim',[-1 1])
% set(hAinterp(2),'ylim',[-7 1])
% set(hA(1),'ylim',[-1 1])
% set(hA(2),'ylim',[-7 1])
% 
% 
% 
% [hx,hy] = format_ticks_v2(hAinterp(1),[],{' '},[],...
%     0,[],[],[],0.02);
% set(hy,'fontsize',font2)
% delete(hx)
% [hx,hy] = format_ticks_v2(hA(1),[],{'$-1$','$-0.5$','$0$','$0.5$','$1$'},[],...
%     -1:0.5:1,[],[],[],0.02);
% set(hy,'fontsize',font2)
% delete(hx)
% 
% 
% [hx,hy] = format_ticks_v2(hAinterp(2),[],{' '},[],...
%     0,[],[],[],0.02);
% set(hy,'fontsize',font2)
% delete(hx)
% % [hx,hy] = format_ticks_v2(hA(2),[],{'$-6$','$-4$','$-2$','$0$'},[],...
% %     -6:2:0,[],[],[],0.0);
% % set(hy,'fontsize',font2)
% % delete(hx)
% 
% %
% set(hA(1),'ycolor',vc);
% set(hA(2),'ycolor',uc);
% 
% 
% ylabh = get(hAinterp(1),'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.05 0 0])
% 
% [hx,hy] = format_ticks_v2(hAinterp(2),{' '},[],0,[],[],[],[],0.005);
% [hx,hy] = format_ticks_v2(hA(1),{' '},[],0,[],[],[],[],0.005);
% [hx,hy] = format_ticks_v2(hA(2),{' '},[],0,[],[],[],[],0.005);
% 
% [hx,hy] = format_ticks_v2(hAinterp(1),{'$0$','$0.25$','$0.5$','$0.75$','$1$'},[],...
%     0:0.25:1,[],[],[],[],0.005);
% set(hx,'fontsize',font2)
% 
% 
% xlabel('$t$ (s)','interpreter','latex','fontsize',font1);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') - [0 .06 0])
% 
% % legend
% hL = legend('$x_{\mathrm{interp}}$','$v_{\mathrm{interp}}$',...
%     '$u_{\mathrm{interp}}$',['$x_{\mathrm{',p.method,'}}$'],...
%     ['$v_{\mathrm{',p.method,'}}$'],['$u_{\mathrm{',p.method,'}}$']);
% set(hL,'orientation','horizontal','interpreter','latex',...
%     'Position',[0.03,0.93,0.95,0.08],'box','off')
% 
% % set(hL,'PlotBoxAspectRatioMode','manual');
% % set(hL,'PlotBoxAspectRatio',[18 1 1]);
% 
% ylabh = get(hAinterp(2),'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [0.04 0 0])
% set(hAinterp(2),'YTick',-6:2:0,'ycolor',uc)
% 
% 
% annotation('textbox',[0.90 0.745 0.1 0.1],'String','$0$','interpreter','latex',...
%     'color',uc,'fontsize',font2,'edgecolor','none')
% 
% annotation('textbox',[0.90 0.545 0.1 0.1],'String','$-2$','interpreter','latex',...
%     'color',uc,'fontsize',font2,'edgecolor','none')
% 
% annotation('textbox',[0.90 0.345 0.1 0.1],'String','$-4$','interpreter','latex',...
%     'color',uc,'fontsize',font2,'edgecolor','none')
% 
% annotation('textbox',[0.90 0.145 0.1 0.1],'String','$-6$','interpreter','latex',...
%     'color',uc,'fontsize',font2,'edgecolor','none')
% 
% %% --- plot error between actual solution
% figI = figI + 1;
% figure(figI); clf(figI);
% % semilogy([10 11],[1 1],'linewidth',2,'color',xc); hold on
% % semilogy([10 11],[1 1],'linewidth',2,'color',vc); hold on
% % semilogy([10 11],[1 1],'linewidth',2,'color',uc); hold on
% % get actual solution values
% [sMat,uMat] = MCW_solution(tarray,p.l);
% % plot position error
% error = abs(interpX(:,1)-sMat(:,1));
% results.maxx = max(error);
% plot(tarray,log10(error),'color',xc,'linewidth',2); hold on
% % plot velocity error
% error = abs(interpX(:,2)-sMat(:,2));
% results.maxv = max(error);
% plot(tarray,log10(error),'color',vc,'linewidth',2); hold on
% % plot control error
% error = abs(interpU-uMat);
% results.maxu = max(error);
% plot(tarray,log10(error),'color',uc,'linewidth',2); hold on
% 
% % label and other plotting stuff
% ylabel('absolute error','interpreter','latex','fontsize',font1);
% 
% set(hL,'interpreter','latex','fontsize',font1-6);
% xlim([p.Tarray(1) p.Tarray(end)]); ylim([-18 0]);
% 
% set(gca,'YTickLabel',[],'YTick',[-15,-10,-5,0])
% 
% [hx,hy] = format_ticks_v2(gca,[],{'$10^{-15}$','$10^{-10}$','$10^{-5}$','$10^{0}$'},[],...
%     [-15,-10,-5,0],[],[],[],0.01);
% set(hy,'fontsize',font2)
% delete(hx)
% 
% [hx,hy] = format_ticks_v2(gca,{'$0$','$0.25$','$0.5$','$0.75$','$1$'},[],...
%     0:0.25:1,[],[],[],[],0);
% set(hx,'fontsize',font2)
% 
% 
% xlabel('$t$ (s)','interpreter','latex','fontsize',font1);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') - [0 .55 0])
% 
% 
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [0.09 0 0])
% 
% % legend
% hL = legend('error $x$','error $v$','error $u$');
% set(hL,'orientation','horizontal','interpreter','latex',...
%     'Position',[0.03,0.93,0.95,0.08],'box','off')
% 
% % set(hL,'PlotBoxAspectRatioMode','manual');
% % set(hL,'PlotBoxAspectRatio',[10 1 1]);
% 
% % results 
% results.Nt = length(U);
% results.NI = length(p.Narray);
% results.f = f;
% results.ferror = abs(f - 4/(9*p.l));
% 
% format shortE
% disp(results)
% 
% % % path = msavename(mfilename('fullpath'),'Saved_Data');
% % % save([path,'BDsol-',p.method,'-Nt_',num2str(results.Nt),'-NI_',num2str(results.NI),...
% % %     '-results','.mat'],'results')
% % 
% % % save2pdf([path,'BDsol-',p.method,'-Nt_',num2str(results.Nt),'-NI_',num2str(results.NI),'-sol.pdf'],1,600)
% % % save2pdf([path,'BDsol-',p.method,'-Nt_',num2str(results.Nt),'-NI_',num2str(results.NI),'-error.pdf'],2,600)
% 
 end