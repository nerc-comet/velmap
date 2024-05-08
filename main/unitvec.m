function [uvec]=unitvec(los,azi,fmt)
%================================================
%function [uvec]=unitvec(los,azi,fmt)
% 
% calculate unit vectors from los and azimuth
% 
% Input:
%  los: line-of-sight angle (positive from vertical)
%  azi: azimuth of the satellite heading vector (positive clockwise from the north)
%  fmt: format (1: roi-pac (default), 2: gamma)
%
% Output:
%  uvec: unit vectors (e,n,u)
%
% see Wang et al., 2007, EPS for the definition
% Hua Wang @ Leeds, 28/10/2009
%================================================

if nargin<3
  fmt=1;
end

n=numel(los);
vlos=reshape(los',n,1);
vazi=reshape(azi',n,1);

%convert from degree to radias
if fmt==1
  vlos = vlos*pi/180.0;
  vazi = -pi/2+vazi*pi/180.0;
else
  vlos = pi/2-vlos;
  vazi = -azi-pi;
end

%unit vectors
uvec(:,1) = [cos(vazi).*sin(vlos)];   %e
uvec(:,2) = [-sin(vazi).*sin(vlos)];  %n
uvec(:,3) = [-cos(vlos)];             %u
