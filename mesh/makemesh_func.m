function [trim]= makemesh_func(x0,x1,y0,y1,dx,dy,dr)
%-------------------------------------------------
%function [trim]= makemesh(x0,x1,y0,y1,dx,dy)
%
% Function to make a regular mesh
% Inputs: 
%   x0/x1/y0/y1: range of x and y axes
%   delx/y: step of x and y axes
%   dr: Maximum fractions of a degree added or subtracted from locations to randomize
% Outputs:
%   trim: structure of the triangular mesh
%     (x/y/trim - x/y coordinates, and connections)
%
% example: trim=makemesh(-119,-116,32,37,0.2,0.2,0.1);
% 
% Hua Wang, 29/09/2015
% Disable saved, added randomisation - Andrew Watson, 19/07/2021
%-------------------------------------------------

%% generate mesh

% generate uniform grids
[y,x]=meshgrid(y0:dy:y1,x0:dx:x1);

% add random shift scaled from 0 to dr
ry = dr .* rand(size(y));
rx = dr .* rand(size(x));
x = x + rx; y = y + ry;

% calculate delaunary triangles
trim.tri = delaunay(x,y);
trim.y=reshape(y,[],1);
trim.x=reshape(x,[],1);

% save trim
