%--------------------------------------------------------------------------
% Copyright Yajur Kumar, 2017. All rights reserved.
% Date : March 13, 2017
%--------------------------------------------------------------------------
% main.m
% Main Program for Calculating Optimal Low-Thrust Motion using
% Legendre-Gauss-Lobatto Pseudospectral Method
% -Main script for solving optimal control problems using 
% multiple-interval pseudospectral methods
% -Define the mesh and method, then solves, then plots solution
% -Note that case 0 can be changed by the user
% -Other cases are based on the technical report of the same name
%--------------------------------------------------------------------------

clc;
clear all;
close all;
study = 4; % 0 indicates user defined mesh

global Re

Re = 6378.14;

% pseudospectral method
p.method = 'LGL'; % either LGL or CGL


% J2 Value
J2    = 1.083e-3;

% Gravitational Parameter
mu    = 398600.4;

global r_0 n_0
% Initial Radius of Reference Orbit
r_0     = 6778;

% Mean Motion
n_0     = sqrt(mu/(r_0^3));

% Maximum Control Input
global Umax Rmin x_final

Umax = 8e-4;
Rmin = -0.25;
x_final = 0.5;

global tf

% Starting and Ending Time
t0 = 0;
tf = 1*pi/n_0;

% switch statement for case studies
switch study
    case 0
        %--- user defined mesh ---%
        p.Tarray = linspace(t0,tf,3); % segment boundaries
        p.Narray = 10*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
        %--- user defined mesh ---%
    case 1
        p.Tarray = [t0,tf]; % segment boundaries
        p.Narray = 6+1; % number of nodes per segment (same size as p.Tarray)
    case 2
        p.Tarray = [t0,tf]; % segment boundaries
        p.Narray = 12+1; % number of nodes per segment (same size as p.Tarray)
    case 3
        p.Tarray = [t0,tf]; % segment boundaries
        p.Narray = 18+1; % number of nodes per segment (same size as p.Tarray)
    case 4
        p.Tarray = linspace(t0,tf,2); % segment boundaries
        p.Narray = 3*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 5
        p.Tarray = linspace(t0,tf,2); % segment boundaries
        p.Narray = 6*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 6
        p.Tarray = linspace(t0,tf,2); % segment boundaries
        p.Narray = 12*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 7
        p.Tarray = linspace(t0,tf,3); % segment boundaries
        p.Narray = 2*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 8
        p.Tarray = linspace(t0,tf,3); % segment boundaries
        p.Narray = 4*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 9
        p.Tarray = linspace(t0,tf,3); % segment boundaries
        p.Narray = 6*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 10
        p.Tarray = linspace(t0,tf,3); % segment boundaries
        p.Narray = 8*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 11
        p.Tarray = linspace(t0,tf,4); % segment boundaries
        p.Narray = 3*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
    case 12
        p.Tarray = linspace(t0,tf,4); % segment boundaries
        p.Narray = 5*ones(1,length(p.Tarray)-1)+1; % number of nodes per segment (same size as p.Tarray)
end

% solve the problem
[t,X,U,f,p]=MCW_solve(p)
MCW_plotConf(t,X,U,f,p)