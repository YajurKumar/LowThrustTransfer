function [C,Ceq] = constraint(s)
C = [];
Ceq = [];

% Equality Constraint: Dynamics
for i=1:s.N
    