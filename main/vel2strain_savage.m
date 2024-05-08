function [strain,eulervec]=vel2strain_savage(trim,vel,vcm,nring,outdir,ndim)
%============================================================
%function [strain,eulervec]=vel2strain_savage(trim,vel,vcm,nring,outdir,ndim)
% forward calculation of strain rate using velocity field in a triangular mesh
% 
% see Savage et al., 2001 JGR, p22,005-22,006 for details
% 
% Input:
%  trim: triangular mesh
%  vel: velocities on the vertices of the mesh
%  vcm: covariance matrix for the velocities
%  nring: ring number
%  outdir: output directory (optional)
%  ndim: dimension of velocities for strain rate estimation (default 2)
%        2 - 2D inversion using A6
%        3 - 3D inversion using A5
%
% Output:
%  strain: strain rates 
%  eulervec: euler vectors
%
% Hua Wang @ Leeds, 22/11/2009
%
% 01/09/2014 HW/CP: 1. add A5 for 3D velocity field
%                   2. add functions to estimate uncertainties
%============================================================

%dimension of the velocity field
%ndim=size(vel,2);
if nargin<6
  ndim=2;
end
if (ndim~=2) && (ndim~=3)
  error('ndim must be 2 or 3');
end

ntri=size(trim.tri,1);
nvtx=length(trim.x);

%initial connectivity matrix
connect = eye(nvtx);
for i=1:ntri
  tri=trim.tri(i,:);
  connect(tri',tri)=1;
end

%connectivity matrix for strain rate calculation
constrain=eye(nvtx);
for ir = 1:nring
  for i=1:nvtx
    constrain(i,:)=constrain(i,:)+sum(connect(constrain(i,:)>0,:),1);
  end
end

%calculate strain rate based on strain rate connectivity matrix
%Savage et al., 2001 JGR, p.22,006, A5
r0 = 6.378e9;  % mean equitorial Earth radius

%2D parameters: [w_theta, w_phi, w_r, e_phiphi, e_thetatheta, e_thetaphi]
%3D parameters: [w_theta, w_phi, w_r, e_phiphi, e_thetatheta, e_thetaphi, U_phi, U_theta, U_r]
m=zeros(nvtx,3*ndim);
phi0=zeros(nvtx,1);
theta0=zeros(nvtx,1);
for i=1:nvtx
  vtxid=find(constrain(i,:)>0);
  
  phi=trim.x(vtxid)*pi/180;        %longitude, convert from degree to rad
  theta=pi/2-trim.y(vtxid)*pi/180; %latitude, change to Colatitude, CW from north
  phi_0(i) = mean(phi);               %centroid longitude
  theta0(i) = mean(theta);           %centroid colatitude
  del_phi = phi-phi_0(i);             %delta longitude
  del_theta = theta-theta0(i);       %delta colatitude

  stheta0=sin(theta0(i));
  ctheta0=cos(theta0(i));

  %design matrix
  n = size(vtxid,2);
  G=zeros(n*ndim,3*ndim);

  if ndim==3 %A5
    G(1:n,3) = del_theta*r0;             %w_r
    G(1:n,4) = stheta0*del_phi*r0;       %e_phiphi 
    G(1:n,6) = del_theta*r0;             %e_thetaphi
    G(1:n,7) = ones(n,1);                %U_phi
    G(1:n,8) = -ctheta0*del_phi;         %U_theta
    G(1:n,9) = -stheta0*del_phi;         %U_r
 
    G(n+1:2*n,3) = -stheta0*del_phi*r0;  %w_r
    G(n+1:2*n,5) = del_theta*r0;         %e_thetatheta
    G(n+1:2*n,6) = stheta0*del_phi*r0;   %e_thetaphi
    G(n+1:2*n,7) = ctheta0*del_phi;      %U_phi
    G(n+1:2*n,8) = ones(n,1);            %U_theta
    G(n+1:2*n,9) = -del_theta;           %U_r
 
    G(2*n+1:3*n,1) = stheta0*del_phi*r0; %w_theta
    G(2*n+1:3*n,2) = -del_theta*r0;      %w_phi
    G(2*n+1:3*n,7) = stheta0*del_phi;    %U_phi
    G(2*n+1:3*n,8) = del_theta;          %U_theta
    G(2*n+1:3*n,9) = ones(n,1);          %U_r

  else %A6
    G(1:n,1) = -ones(n,1);               %w_theta
    G(1:n,2) = -ctheta0*del_phi;         %w_phi
    G(1:n,3) =  del_theta;               %w_r
    G(1:n,4) =  stheta0*del_phi;         %e_phiphi
    G(1:n,6) =  del_theta;               %e_thetaphi
 
    G(n+1:2*n,1) = -ctheta0*del_phi;     %w_theta
    G(n+1:2*n,2) =  ones(n,1);           %w_phi
    G(n+1:2*n,3) = -stheta0*del_phi;     %w_r
    G(n+1:2*n,5) =  del_theta;           %e_thetatheta
    G(n+1:2*n,6) =  stheta0*del_phi;     %e_thetaphi
 
    G=G*r0;
  end

  %velocity
  dvel = reshape(vel(vtxid',1:ndim),ndim*n,1);
  dvel(n+1:2*n) = -dvel(n+1:2*n); %southward for v_n, corresponding to colatitude

  %vcm of velocity
  sel=[vtxid,nvtx+vtxid];
  if ndim==3
    sel=[sel,2*nvtx+vtxid];
  end
  vcm_tmp = full(vcm(sel,sel'));

  %least-squares solution
  [m(i,:),~,mse,covm(i,:,:)]=lscov(G,dvel,vcm_tmp);
  %covm(i,:,:)=covm(i,:,:)/mse;

  %uncertainty for strain and rotation
  %s_w: [s_w_theta,s_w_phi,s_w_r];
  %s_deps: [s_e_phiphi,s_e_thetatheta,s_e_thetaphi];
  s_w(i,1:3)=sqrt(diag(squeeze(covm(i,1:3,1:3))));
  s_deps(i,1:3)=sqrt(diag(squeeze(covm(i,4:6,4:6))));

  % -m(:,6): from colatitude to latitude
  % covm(6,6): without changed
  m(i,6)=-m(i,6);
  covm(i,:,6)=-covm(i,:,6);
  covm(i,6,:)=-covm(i,6,:);
end

%calculate Euler vector, Savage et al., 2001, p22,006, A7
[pole,s_pole]=eulertrans(m(:,1:3),theta0,phi0,covm(:,1:3,1:3));
%final euler vectors (NB: Carolina outputs -m(:,1:3))
eulervec=[trim.x,trim.y,m(:,1:3),pole,s_w,s_pole];

%calculate eps_rr Savage 2001 JGR pg 22,006, assuming vu=0.25
%eps_rr = (-m(:,4)-m(:,5))/3;

%calculate principal strain rates
[peps,s_peps]=prinstrain(m(:,4:6),1,covm(:,4:6,4:6));
%final strain
strain=[trim.x,trim.y,m(:,4:6),peps,s_deps,s_peps];

if nargin>2
  %format: lon,lat,e_phiphi,e_thetatheta,e_thetaphi,dilat,maxshear,eps1,eps2,alpha,dE,...
  %                sigma of the above parameters
  save(strcat(outdir,'strain_savage','_nring',num2str(nring),'.dat'),'strain','-ASCII');
  save(strcat(outdir,'euler_savage','_nring',num2str(nring),'.dat'),'eulervec','-ASCII');
end
