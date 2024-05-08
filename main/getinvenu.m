function [invenu] = getinvenu(gps,insar)
%======================================================
%function [invenu] = getinvenu(gps,insar)
%
% Get velocity inversion parameters from GPS & InSAR dataset
%                                                      
% INPUT:                                               
%   gps: gps dataset
%   insar: insar dataset (optional)
% OUTPUT:                                              
%   invenu: velocity inversion parameters (inv_e, inv_n, inv_u)                         
%                                                      
% Hua Wang, 28/05/2011                         
%
% 16/03/2015 HW: support loading multiple gps files
%======================================================

if nargin<2
  insar=[];
end

%determine invenu by gps data alone
ngf=length(gps);
invenu=[0 0 0];
for i=1:ngf
  invenu=invenu+gps(i).invenu;
end
  
%determine invenu by insar data (e/n/u)
ninsar=length(insar);
for is=1:ninsar
  invenu=invenu+insar(is).proc.invenu;
end
invenu(invenu~=0)=1;
