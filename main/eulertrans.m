function [pole,s_pole]=eulertrans(m,theta0,phi0,vcm)
%============================================================
%function [pole,s_pole]=eulertrans(m,theta0,phi0,vcm)
%calcualte Euler pole from omega_theta, omega_phi, omega_r
%
% Inputs:
%   m:      matrix [omega_theta, omega_phi, omega_r]
%   theta0: centoid colatitude
%   phi0:   centoid longitude
%   vcm:    covariance matrix of m
%
% Output:
%   pole:  Euler pole [phi_p,theta_p,omega]
%   s_pole: 1-sigma of Euler pole
%
%ref: Savage et al., 2001, p22,006, A7
%
% Hua Wang, 1/9/2015
%============================================================

omega_theta=m(:,1);
omega_phi=m(:,2);
omega_r=m(:,3);
omega = sqrt(omega_theta.^2+omega_phi.^2+omega_r.^2); %rotation rate
y = omega_r.*cos(theta0) - omega_theta.*sin(theta0);
theta_p = acos(y./omega)*180/pi; %colatitude of rotation pole
A = omega_r.*sin(theta0).*sin(phi0) + omega_theta.*cos(theta0).*sin(phi0) + omega_phi.*cos(phi0);
B = omega_r.*sin(theta0).*cos(phi0) + omega_theta.*cos(theta0).*cos(phi0) - omega_phi.*sin(phi0);
phi_p=atan2(A,B)*180/pi;         %longitude of rotation pole

theta_p=90-theta_p;          %latitude of rotation pole
id=find(theta_p<0);
theta_p(id)=-theta_p(id);
phi_p(id)=phi_p(id)-180;
omega(id)=-omega(id);        %ccw positive, unit in rad/yr

pole=[phi_p,theta_p,omega];

%-----------
%uncertainty estimates, linearise from A7
%-----------

%linearise of phi
J1(:,1) = (B.*sin(phi0) - A.*cos(phi0)).*cos(theta0)./(A.^2 + B.^2); %derivative of omega_theta
J1(:,2) = (B.*cos(phi0) + A.*sin(phi0))./(A.^2 + B.^2);               %derivative of omega_phi
J1(:,3) = (B.*sin(phi0) - A.*cos(phi0)).*sin(theta0)./(A.^2 + B.^2); %derivative of omega_r

%linearise of theta
J2(:,1)=-sin(theta0) - y.*omega_theta./omega.^2; %derivative of omega_theta
J2(:,2)=-y.*omega_phi./omega.^2;                  %derivative of omega_phi
J2(:,3)= cos(theta0) - y.*omega_r./omega.^2;     %derivative of omega_r
%derivate of arc-cosine
z = -1./sqrt(1 - (y./omega).^2);
J2=J2.*repmat(z./omega,1,3);

%linearise of omega
J3=m./repmat(omega,1,3);

n=size(m,1);
s_pole=zeros(n,3);
for i=1:n
  ivcm=squeeze(vcm(i,:,:));
  J=[J1(i,:);J2(i,:);J3(i,:)];
  p_vcm=J*ivcm*J';
  s_pole(i,:)=sqrt(diag(p_vcm));
end
