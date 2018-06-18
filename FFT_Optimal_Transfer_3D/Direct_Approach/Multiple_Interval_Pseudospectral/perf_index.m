%--------------------------------------------------------------------------
% perf_index.m
% Performance Index for the Optimization Problem
% Calculates the Lagrange term for the performance index
%--------------------------------------------------------------------------
% PI = perf_index(x,p)
% x          : optimization states vector
% p          : parameter structure
% PI         : Performace Index Value
% h          : t(N_t)-t(0)
% w_CGL(i)   : weight at the ith step as per CGL
% w_LGL(i)   : weight at the ith step as per LGL
% Lag_Energy : Lagrange term for minimum energy
%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%
%--------------------------------------------------------------------------

function PI = perf_index(x,p)
PI  =  0; % Initializing PI

% Evaluating the value inside the integral

if p.Choice == 1
    for i=0:N_t
    PI  =  PI + (h/2)*w_CGL(i)*Lag_Energy(i);
end
end
else if p.Choice == 2
        for i=0:N_t
        PI  =  PI + (h/2)*w_LGL(i)*Lag_Energy(i);
        end
    end
    