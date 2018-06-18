%--------------------------------------------------------------------------
% scaling.m
% scaling function that changes variable bounds to be between 0 and 1 and
% vice versa
% needs vector inputs
%--------------------------------------------------------------------------
% x_out = scaling(x,l,u,type)
%    x: input variable vector
%    l: unscaled lower bound
%    u: unscaled upper bound
% type: scale (1) or unscale (2)
% x_out: output variable vector
%--------------------------------------------------------------------------
% Author: Daniel R. Herber, Graduate Student, University of Illinois at
% Urbana-Champaign
% Date: 06/04/2015
%--------------------------------------------------------------------------
function x_out = scaling(x,l,u,type)
    if type == 1
        x_out = (x-l)./(u-l); % scale
    elseif type == 2
        x_out = l + x.*(u-l); % 
    else
        error('wrong type')
    end
end