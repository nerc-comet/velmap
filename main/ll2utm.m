function xy=ll2utm(ll,origin)
%----------------------------------------
%function xy=ll2utm(ll,origin)
%  Converts from longitude and latitude to utm
% Inputs:
%   ll: (lon,lat) - n by 2 matrix
%   origin: (lon,lat) of origin
% Outputs:
%   xy: (xe,yn) w.r.t the origin - n by 2 matrix
%
% Hua Wang, 17 Oct 2013
%
% Note: WGS84 ellipsoid is used
%----------------------------------------

%Set ellipsoid constants (WGS84)
a=6378137.0;
e2=0.00673949674227;
e4=e2^2;
e6=e2^3;
e8=e2^4;
c=a*sqrt(1+e2);
b0=c*(1-3*e2/4+45*e4/64-175*e6/256+11025*e8/16384);
b1=b0-c;
b2=c*(15*e4/32-175*e6/384+3675*e8/8192);
b3=c*(-35*e6/96+735*e8/2048);
b4=c*315*e8/1024;

%Convert to radians
ll=[ll;origin]; %append origin to ll
ll(:,1)=ll(:,1)-origin(1);
ll=double(ll)*pi/180;
cos1=cos(ll(:,2));
sin1=sin(ll(:,2));

%Formula 7-107, P18
y=b0*ll(:,2)+(b1+b2*cos1.^2+b3*cos1.^4+b4*cos1.^6).*cos1.*sin1;

%Formula 8-97, P82
eta2=e2*cos1.^2;
N=c./sqrt(1+eta2);
m=cos1.*ll(:,1);
t=tan(ll(:,2));
t2=t.^2;
xy(:,1)=N.*(m+(1-t2+eta2).*(m.^3)/6+(5-18*t2+t2.^2+14*eta2-58*eta2.*t2)/120.*(m.^5));
xy(:,2)=y+N.*t.*(0.5*(m.^2)+(5-t2+9*eta2+4*eta2.^2).*(m.^4)/24+(61-58*t2+t2.^2)/720.*(m.^6));

%Convert to km
xy=xy/1000;

%remove origin from xy
npt=size(xy,1);
xy(:,2)=xy(:,2)-xy(npt,2);
xy(npt,:)=[];
