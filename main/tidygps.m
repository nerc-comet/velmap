function [gps]=tidygps(trim,gps)
%==========================================
%function [gps]=tidygps(trim,gps)
%
% Remove some gps data which are outside the mesh
%
% Input:
%   trim: triangular mesh
%   gps:  gps data
%
% Output:
%   gps:  tidyed gps data
%
% Hua Wang @ Uni Leeds, 24/09/2009
%
% 19/08/2011 HW: change the structure of GPS for time series
% 27/03/2010 HW: support 3D velocities
%=========================================

%coarse tidy up gps
xr=minmax(trim.x');
yr=minmax(trim.y');
lon=[gps.lon];  %lon/lat of all gps stations
lat=[gps.lat];  %lon/lat of all gps stations
tidy=find(lon<xr(1) | lon>xr(2) | lat<yr(1) | lat>yr(2));
gps(tidy)=[];

%precise tidy up gps
ntri=length(trim.tri);    %triangular number
ngps=length(gps);      
%design matrix dimension
locate=zeros(ngps,1);

for itri = 1:ntri
  %vert coordinates for each triangular
  x=trim.x(trim.tri(itri,:)); %x coordinates of the three vertex
  y=trim.y(trim.tri(itri,:)); %y coordinates of the three vertex

  for is=1:ngps
    if locate(is)==0
      pos=intri([gps(is).lon,gps(is).lat],[x y]);
      if pos~=0
        locate(is)=1;
      end
    end
  end
end
tidy=find(locate==0);
gps(tidy)=[];
