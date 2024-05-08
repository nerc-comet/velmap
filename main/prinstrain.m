function [peps,s_peps]=prinstrain(eps,prin,vcm)
%============================================================
%function [peps,s_peps]=prinstrain(eps,prin,vcm)
% forward calculation of principal strain rate
% 
% Input:
%  eps: 2-d strain rate tensor (epsilon_xx,epsilon_yy,epsilon_xy)
%  prin: method to calculate principal strain rate (default 1 eig)
%  vcm: covariance matrix of eps
%
% Output:
%  peps: principal strain rates
%  s_peps: 1-sigma uncertainty of principal strain rates
%
% Hua Wang @ Leeds, 19/11/2009
%
% 02/09/2015 HW: incorporate two methods for prin strain
% 01/09/2015 HW/CP: uncertainty estimates
%============================================================

if nargin<2
  prin=1;
end

n=size(eps,1);

gamma1=0.5*(eps(:,1)-eps(:,2));  %(exx-eyy)/2
gamma2=eps(:,3);                 %exy

%dilatation and maximum shear strain rate
maxshear = sqrt(gamma1.*gamma1+gamma2.*gamma2); %p118, Turcott&Schubert v3.
dilat=eps(:,1)+eps(:,2);         %exx+eyy, p112, Turcott&Schubert v3.

%principal strain rate and direction
if prin==2
  %max/min principal strain rate (eq 2.129, p118, Turcott&Schubert v3.)
  eps1 = 0.5*dilat+maxshear;
  eps2 = 0.5*dilat-maxshear;
  alpha = -atan2(gamma2,gamma1)*90/pi; %clockwise from north, not correct here?
  alpha(alpha<0)=alpha(alpha<0)+180; % 0-180
else
  eps1=zeros(n,1);
  eps2=zeros(n,1);
  alpha=zeros(n,1);
  for i=1:n
    A=[eps(i,1) eps(i,3); eps(i,3) eps(i,2)];  %strain rate tensor
    [V,D]=eig(A);                              %eignvectors and eignvalues
    [eps1(i),imax]=max(diag(D));               %e1: maximum eignvalue (mostly extension)
    [eps2(i),imin]=min(diag(D));               %e2: minimum eignvalue (mostly compression)
    vmin=V(:,imin);                            %using e2 to determine principal axis, consistent with GMT
    alpha(i)=90-atan2(vmin(2),vmin(1))*180/pi; %clockwise from north
  end
  %alpha(alpha<0)=alpha(alpha<0)+360;    %to 0-360
  %alpha(alpha>180)=alpha(alpha>180)-180;  %to 0-180
end

%second invariant of the strain rate tensor, revised on 05/10/2010
%England and McKenzie, 1982, formular (10), without divided by 2
dE=sqrt(eps(:,1).^2 + eps(:,2).^2 + eps(:,3).^2+eps(:,1).*eps(:,2));
peps=[dilat,maxshear,eps1,eps2,alpha,dE];

%---------
%added on 1/9/2015 for uncertainty estimates
%---------
if nargin>2
  s_peps=zeros(n,6); % s_dilat,s_maxshear,s_eps1,s_eps2,s_alpha,s_dE
  J=zeros(6,3);
  for i=1:n
    %sigma of dilatation, exx+eyy=e1+e2
    J(1,:)=[1 1 0];      %build Jacobian matrix, J
 
    %sigma of maximum shear, derived from eq. 2.129, p. 119, Turcott&Schubert v3.
    J(2,:)=[gamma1(i)/2, -gamma1(i)/2, gamma2(i)]/maxshear(i);   %build Jacobian matrix, J
 
    %sigma of eps1/eps2/alpha, derived from eq. 2.129, p. 119, Turcott&Schubert v3.
    J(3,:)= J(2,:)+[0.5 0.5 0];
    J(4,:)=-J(2,:)+[0.5 0.5 0];

    %sigma of alpha, derived from eq. 2.128, p. 119, Turcott&Schubert v3.
    %here is different from Hammond, his is: J(5,:)=[eyy, -eyy, -2*gamma1(i)];
    J(5,:)=[gamma2(i), -gamma2(i), -2*gamma1(i)];
    J(5,:)=J(5,:)/(maxshear(i)^2)*45/pi; %alpha convert to degree
 
    %sigma of dE
    J(6,:)=[eps(i,1)+0.5*eps(i,2), 0.5*eps(i,1)+eps(i,2), eps(i,3)]/dE(i);
 
    %vcm of all the strains
    vcm_peps=J*squeeze(vcm(i,:,:))*J';
    s_peps(i,:)=sqrt(diag(vcm_peps));
  end
end
