%--------------------------------------------------------------------------
% Method_Pseudospectral.m
% Attempt to solve the Bryson-Denham problem using a pseudospectral method
% (namely LGL nodes and Gaussian quadrature) 
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% https://github.com/danielrherber/optimal-control-direct-method-examples
%--------------------------------------------------------------------------
function Method_Pseudospectral
    
	mu    = 398600.4;
	R     = 6778.14;
	
	global N;
	N = sqrt(mu/(R^3));
	
	% problem parameters
    p.ns = 4; p.nu = 2; % number of states and controls
    p.t0 = 0; p.tf = 1; % time horizon
	
	
	
	global x0 Vx0 y0 Vy0 xf Vxf yf Vyf

% Initial Conditions
p.x0  = 10;   % in km
p.Vx0 = 0.01; % in km/s
p.y0  = 2;   % in km
p.Vy0 = 0;    % in km/s

% Final Conditions
p.xf  = 0;    % in km
p.Vxf = 0;    % in km/s
p.yf  = 0;    % in km
p.Vyf = 0;    % in km/s
	
    % p.y10 = 0; 
	% p.y1f = 0; 
	% p.y20 = 1; 
	% p.y2f = -1; % boundary conditions
    
	% p.l = 1/9;
    
	% direct transcription parameters
    p.nt = 10; % number of node points
%     p.nt = 50; % number of node points
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
    options = optimoptions(@fmincon,'display','iter','MaxFunctionEvaluations',1e5); % options
    
	% solve the problem
    x = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],@(x) constraints(x,p),options);
    
	% obtain the optimal solution
    X   = x(p.X); 
	VX  = x(p.VX);
    Y   = x(p.Y);
    VY	= x(p.VY)
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	
	p.t = (p.tau*(p.tf-p.t0) + (p.tf+p.t0))/2; % unscale time horizon
    % plots
    Plots(X,VX,Y,VY,UX,UY,p,'Pseudospectral')
end

% objective function
% function f = objective(x,p)
    % u = x(p.ui); % extract
    % L = u.^2; % integrand
    % f = (p.tf-p.t0)/2*dot(p.w,L)/2; % calculate objective
% end

function fobj = objective(x,p)
    fobj = 0; % initialize
    % determine integral value
	
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	L = UX.^2 + UY.^2; % integrand
    fobj = (p.tf-p.t0)/2*dot(p.w,L)/2; % Gaussian quadrature  
end


% constraint function
% function [c,ceq] = constraints(x,p)

    % y1 = x(p.y1i); 
	% y2 = x(p.y2i); 
	% u = x(p.ui); % extract
	
    % Y = [y1,y2]; F = (p.tf-p.t0)/2*[y2,u]; % create matrices (p.nt x p.ns)
    
	% ceq1 = y1(1) - p.y10; % initial state conditions
    % ceq2 = y2(1) - p.y20;
    % ceq3 = y1(end) - p.y1f; % final state conditions
    % ceq4 = y2(end) - p.y2f;
    % ceq5 = p.D*Y - F; % defect constraints matrix form (pseudospectral)
    
	% c1 = y1 - p.l; % path constraints
    % c = c1; 
	% ceq = [ceq1;ceq2;ceq3;ceq4;ceq5(:)]; % combine constraints
% end


function [c,ceq] = constraints(x,p)

    % initialize constraints
    c = []; ceq = []; % contin = []; in_contin = [];
    % extract from optimization vector
   
	X   = x(p.X); 
	VX  = x(p.VX);
    Y   = x(p.Y);
    VY	= x(p.VY)
	UX  = x(p.UX); % extract
    UY  = x(p.UY);
	X = [X,VX,Y,VY];
	U = [UX,UY];
	
	global N;
	
    %--- equality constraints
    
	F = [x(2), 3*N*N*x(1)+2*N*x(4)-x(6), x(4), -2*N*x(2)-x(8), -3*N*N*x(6), -x(5)+2*N*x(8), 0, -2*N*x(6)-x(7)];
	cons = p.D*X - F; % defect constraints matrix form (pseudospectral)
	
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
    
    % Rx = X(:,1);
    % Ry = X(:,3);
    % Rz = X(:,5);
    
    % Ux = U(:,1);
    % Uy = U(:,2);
    % Uz = U(:,3);
    
    
    % Rnorm = (Rx.^2 + Ry.^2).^(0.5);
    % Unorm = (Ux.^2 + Uy.^2).^(0.5);
    % c = [c;Unorm-Umax];
    % c = [c;-Rnorm+Rmin];
    

     
end