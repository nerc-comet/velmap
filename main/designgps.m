function [interpmat]=designgps(trim,gps,invenu)
%============================================================
%function [interpmat]=designgps(trim,gps,invenu)
% 
% Design interpolation matrix for gps data
% 
% Input:
%  trim: triangular mesh structure (tri,x,y)
%  gps:  gps data structure
%  invenu: inverse parameters
%
% Output:
%  interpmat: interpolation design matrix
%
% Hua Wang @ Leeds, 23/09/2009
%
% 02/09/2011 HW: remove uninversed components
% 27/03/2010 HW: support 3D velocities
%============================================================
if nargin<3
  invenu=[1 1 1];
end

ntri=length(trim.tri);    %triangular number
nvtx=length(trim.x);      %vertex number
ngps=length(gps);         %gps stations number

%design matrix dimension
kernel=sparse(ngps,nvtx);
locate=zeros(ngps,1,'int8');

for itri = 1:ntri
  
  %vertex coordinates for each triangular
  x=trim.x(trim.tri(itri,:)); %x coordinates of the three vertex
  y=trim.y(trim.tri(itri,:)); %y coordinates of the three vertex

  %improve the efficiency by rectangular test first
  xr=minmax(x');
  yr=minmax(y');
  lon=[gps.lon];
  lat=[gps.lat];
  ptin=find(lon>=xr(1) & lon<=xr(2) & lat>=yr(1) & lat<=yr(2));
  nleft=length(ptin);
  for i=1:nleft
    is=ptin(i);
    if locate(is)==0
      pt=[gps(is).lon, gps(is).lat];
      tri=[x,y];
      pos=intri(pt,tri);
      if pos~=0
        %column numbers in the design matrix
        icol=trim.tri(itri,:);
        %interpolation function (England & Molnar, 2005, P5, Fm5-8)
        kernel(is,icol)=interpk(pt,tri);
      
        %so that this pixel will not be evaluated for another triangular
        locate(is)=1;
      end
    end
  end
end

if nnz(locate)<ngps
  error('not all points have been evaluated for triangular search')
end

%make the final design matrix
%could be a bug if GPS components < invenu ????  HW
ncomp=sum(invenu);
interpmat=kernel;
for i=2:ncomp
  interpmat=blkdiag(interpmat,kernel);
end
