%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% MCW_solve.m
%--------------------------------------------------------------------------
% [t,X,U,f,p] = PS_solve(p)
% p: parameter structure
% t: time vector
% X: state matrix
% U: control matrix
% f: objective function value
% States:
% X1 : x
% X2 : x'
% X3 : y
% X4 : y'
% X5 : z
% X6 : z'
%--------------------------------------------------------------------------

function [t,X,U,f,p] = MCW_solve_loop(p)
global Umax Rmin x_final
%--- problem specific parameters ---%
p.nt = sum(p.Narray); % total number of nodes
p.ns = 6; % number of states
p.nu = 3; % number of controls

% Boundary Conditions
p.X1_0 = 18; 
p.X2_0 = 0.001; 
p.X3_0 = 12; 
p.X4_0 = 0.003; 
p.X5_0 = 10; 
p.X6_0 = 0.001; 
p.X1_f = x_final; 
p.X2_f = 0; 
p.X3_f = x_final; 
p.X4_f = 0; 
p.X5_f = x_final; 
p.X6_f = 0; 

% variable bounds
%XL = -Inf*ones(p.nt,p.ns); % state lower
%XU = Inf*ones(p.nt,p.ns); % state upper
%UL = -Umax*ones(p.nt,p.nu); % control lower
%UU = Umax*ones(p.nt,p.nu); % control upper
%p.xL = [reshape(XL,[],1);reshape(UL,[],1)]; % combine to column vector
%p.xU = [reshape(XU,[],1);reshape(UU,[],1)]; % combine to column vector

p.deriv = 'MCW_deriv'; % derivative function name
 %p.deriv = 'CW_deriv'; % derivative function name
% p.deriv = 'CW_deriv'; % derivative function name
p.lagrange = 'MCW_lagrange'; % Lagrange term function name
p.mayer = 'MCW_mayer'; % Mayer term function name
%p.path = 'MCW_path'; % path constraint function name
p.boundary = 'MCW_boundary'; % boundary constraint function name
p.U_max = 'U_max'; % constraint on the applied control input

%--- mesh parameters ---%
p.cumN = [0,cumsum(p.Narray)]; % segment ending indices

% build segment information
for i = 1:length(p.Narray)
    if strcmp(p.method,'LGL')
        p.tau{i,1} = LGL_nodes(p.Narray(i)-1); % calculate scaled node locations
        p.w{i,1} = w_LGL(p.tau{i}); % segment Gaussian quadrature weights
        p.D{i,1} = LGL_Dmatrix(p.tau{i}); % segment differentiation matrix
    elseif strcmp(p.method,'CGL')
        p.tau{i,1} = CGL_nodes(p.Narray(i)-1); % calculate scaled node locations
        p.w{i,1} = w_CGL(p.tau{i}); % segment Gaussian quadrature weights
        p.D{i,1} = CGL_Dmatrix(p.tau{i}); % segment differentiation matrix
    end
    p.t0(i,1) = p.Tarray(i); % segment initial time
    p.tf(i,1) = p.Tarray(i+1); % segment final time
    p.t{i,1} = ((1-p.tau{i})*p.t0(i) + (1+p.tau{i})*p.tf(i))/2; % segment unscaled time values
end

% initial guess
% linear guess functions
fX1_0 = @(t) (p.X1_f - p.X1_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X1_0;
fX2_0 = @(t) (p.X2_f - p.X2_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X2_0;
fX3_0 = @(t) (p.X3_f - p.X3_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X3_0;
fX4_0 = @(t) (p.X4_f - p.X4_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X4_0;
fX5_0 = @(t) (p.X5_f - p.X5_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X5_0;
fX6_0 = @(t) (p.X6_f - p.X6_0)/(p.t{end}(end) - p.t{1}(1))*t + p.X6_0;

fU1_0 = @(t,u10,u1f) (u1f - u10)/(p.t{end}(end) - p.t{1}(1))*t + u10;
fU2_0 = @(t,u20,u2f) (u2f - u20)/(p.t{end}(end) - p.t{1}(1))*t + u20;
fU3_0 = @(t,u30,u3f) (u3f - u30)/(p.t{end}(end) - p.t{1}(1))*t + u30;

X10 = []; 
X20 = []; 
X30 = []; 
X40 = []; 
X50 = []; 
X60 = []; 

U10 = [];
U20 = [];
U30 = [];


for i = 1:length(p.Narray)
    X10 = [X10;fX1_0(p.t{i})]; % X1 initial guess
    X20 = [X20;fX2_0(p.t{i})]; % X2 initial guess
    X30 = [X30;fX3_0(p.t{i})]; % X3 initial guess
    X40 = [X40;fX4_0(p.t{i})]; % X4 initial guess
    X50 = [X50;fX5_0(p.t{i})]; % X5 initial guess
    X60 = [X60;fX6_0(p.t{i})]; % X6 initial guess
    
    U10 = [U10;fU1_0(p.t{i},0,0)]; % X control initial guess
    U20 = [U20;fU2_0(p.t{i},0,0)]; % Y control initial guess
    U30 = [U30;fU3_0(p.t{i},0,0)]; % Z control initial guess
end
x0 = [reshape(X10,[],1);reshape(X20,[],1);reshape(X30,[],1);reshape(X40,[],1);reshape(X50,[],1);reshape(X60,[],1);reshape(U10,[],1);reshape(U20,[],1);reshape(U30,[],1)]; % combine to column vector

%--- optimization ---%
% optimization routine options
 options = optimoptions(@fmincon,'Display','Iter','Algorithm','sqp',...
    'MaxFunctionEvaluations',Inf,'MaxIterations',10000,'FunctionTolerance',1e-8,'StepTolerance',1e-8,'OptimalityTolerance',1e-8); % interior-point

%options = optimoptions(@fmincon,'display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',10000,'StepTolerance',1e-15); % options
tic % start timer

    %[x,f] = fmincon(@(x) objective(x,p),x0,[],[],[],[],p.xL,p.xU,...
     %   @(x) nonlincon(x,p),options);
    
    [x,f,exitflag,output] = fmincon(@(x) objective(x,p),x0,[],[],[],[],[],[],...
        @(x) nonlincon(x,p),options);
toc % end timer
p.iterate  = output.iterations;
p.funcount = output.funcCount;
p.constvio = output.constrviolation;
p.stepsz   = output.stepsize;
p.algo     = output.algorithm;
p.J_cost = f;
p.msg      = output.message;
p.time     = toc;
% extract from optimization vector
    X = reshape(x(1:p.ns*p.nt),p.nt,[]); % states
    U = reshape(x(p.ns*p.nt+1:end),p.nt,[]); % controls
    % Starting and Ending Markers
    
    xstart = 1;
    xend   = xstart + (p.ns/6)*p.nt - 1;
    vxstart = xend + 1;
    vxend   = vxstart + (p.ns/6)*p.nt - 1;
    ystart = vxend + 1;
    yend   = ystart + (p.ns/6)*p.nt - 1;

    vystart= yend + 1;
    vyend  = vystart + (p.ns/6)*p.nt - 1;
    zstart= vyend + 1;
    zend  = zstart + (p.ns/6)*p.nt - 1;
    vzstart= zend + 1;
    vzend  = vzstart + (p.ns/6)*p.nt - 1;
    
    uxstart = p.ns*p.nt+1;
    uxend   = uxstart + ((p.nu/3)*p.nt) - 1;
    uystart = uxend + 1;
    uyend   = uystart + ((p.nu/3)*p.nt) - 1;
    uzstart = uyend + 1;
    uzend   = uzstart + ((p.nu/3)*p.nt) - 1; 
    
    Rx = reshape(x(xstart:xend),p.nt,[]);
    Vx = reshape(x(vxstart:vxend),p.nt,[]);
    Ry = reshape(x(ystart:yend),p.nt,[]);
    
    Vy = reshape(x(vystart:vyend),p.nt,[]);
    Rz = reshape(x(zstart:zend),p.nt,[]);
    Vz = reshape(x(vzstart:vzend),p.nt,[]);
    
    Ux = reshape(x(uxstart:uxend),p.nt,[]);
    Uy = reshape(x(uystart:uyend),p.nt,[]);
    Uz = reshape(x(uzstart:uzend),p.nt,[]);

% unscaled time vector
t = cell2mat(p.t);

end