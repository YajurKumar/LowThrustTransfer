
function Method_Pseudospectral
    
clear all;
close all;
clc;
	mu    = 398600.4;
	R     = 6778.14;
	
	global N;
	N = sqrt(mu/(R^3));
	
	% problem parameters
    p.ns = 4; p.nu = 2; % number of states and controls
    p.t0 = 0; p.tf = 1000; % time horizon
	
	
	
	global x0 Vx0 y0 Vy0 xf Vxf yf Vyf

% Initial Conditions
p.x0  = 100;   % in km
p.Vx0 = 0.01; % in km/s
p.y0  = 200;   % in km
p.Vy0 = 0;    % in km/s

% Final Conditions
p.xf  = 0;    % in km
p.Vxf = 0;    % in km/s
p.yf  = 0;    % in km
p.Vyf = 0;    % in km/s

    
	% direct transcription parameters
    p.nt = 50; % number of node points
    p.tau = LGL_nodes(p.nt-1); % scaled time horizon
    p.D =  LGL_Dmatrix(p.tau); % for defect constraints
    p.w = w_LGL(p.tau); % for gaussian quadrature
	
    % discretized variable indices in x = [x,Vx,y,Vy,Ux,Uy];
    p.X   = 1:p.nt; 
	p.VX  = p.nt+1:2*p.nt; 
	p.Y   = 2*p.nt+1:3*p.nt;
	p.VY  = 3*p.nt+1:4*p.nt;
	p.UX  = 4*p.nt+1:5*p.nt;
	p.UY  = 5*p.nt+1:6*p.nt;
	
	
    x0 = zeros(p.nt*(p.ns+p.nu),1); % initial guess (all zeros)
    options = optimoptions(@fmincon,'display','iter','MaxFunctionEvaluations',1e6); % options
    
	% solve the problem
    x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
    
	% obtain the optimal solution
    X   = x(p.X); 
	VX  = x(p.VX);
    Y   = x(p.Y);
    VY	= x(p.VY);
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	
	p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon
   figure
   subplot(3,1,1)
   plot(p.t,X,p.t,Y);
   xlabel('time (s)');
   ylabel('distances (km)');
   legend('x','y')
   title('Relative Distance');
   U = sqrt(UX.^2+UY.^2);
   %figure
   subplot(3,1,2)
   plot(p.t,U);
   xlabel('time (s)');
   ylabel('Specific Thrust (km/s^2)')
   title('Control Input')
   subplot(3,1,3)
   plot(X,Y);
   xlabel('x (km)');
   ylabel('y (km)');
   title('Geometry of the orbit')
end


function fobj = objective(x,p)
    
    % determine integral value
	
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	L = UX.^2 + UY.^2; % integrand
    fobj = (p.tf-p.t0)/2*dot(p.w,L)/2; % Gaussian quadrature  
end




function [c,ceq] = constraints(x,p)

    % initialize constraints
    c = []; ceq = []; % contin = []; in_contin = [];
    % extract from optimization vector
   
	X   = x(p.X); 
	VX  = x(p.VX);
    Y   = x(p.Y);
    VY	= x(p.VY);
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	X_vec = [X,VX,Y,VY];
	U = [UX,UY];
	
	global N;
	
    %--- equality constraints
    
	F = (p.tf-p.t0)/2*[VX, 3*N*N*X+2*N*VY+UX, VY, -2*N*VX+UY];
    cons = p.D*X_vec - F; % defect constraints matrix form (pseudospectral)
	
	ceq1 = X(1) - p.x0;
	ceq2 = VX(1) - p.Vx0;
	ceq3 = Y(1) - p.y0;
	ceq4 = VY(1) - p.Vy0;
    ceq5 = X(end) - p.xf;
    ceq6 = VX(end) - p.Vxf;
    ceq7 = Y(end) - p.yf;
    ceq8 = VY(end) - p.Vyf;
	
	ceq = [ceq1; ceq2; ceq3; ceq4; ceq5; ceq6; ceq7; ceq8; cons(:)];
	
	% % determine defect constraints for each segment
    % for i = 1:length(p.Narray)
        % Xi = X(p.cumN(i)+1:p.cumN(i+1),:); % states for segment i
        % Ui = U(p.cumN(i)+1:p.cumN(i+1),:); % controls for segment i
        % D = p.D{i}; % differentiation matrix
        % F = Fmatrix(p.tau{i},Xi,Ui,p.t0(i),p.tf(i),p); % differential matrix
        % defect{i,1} = D*Xi - F; % defect constraints
    % end
    % ceq = [ceq;reshape(cell2mat(defect),[],1)];
    
    % boundary constraints
    % boundary = eval([p.boundary,'(X,U,p)']);
    % ceq = [ceq;reshape(boundary,[],1)];
    
   
    
    Umax = 1e-6;
    Rmin = 0;
    Rnorm = (X.^2 + Y.^2).^(0.5);
    Unorm = (UX.^2 + UX.^2).^(0.5);
    c = [c;Unorm-Umax];
    c = [c;-Rnorm+Rmin];
    

     
end