function [strain]=vel2strain_tri(trim,vel,outdir)
%============================================================
%function [strain]=vel2strain_tri(trim,vel,outdir)
% forward calculation of strain rate using velocity field in
% a triangular mesh
% 
% Input:
%  trim: triangular mesh
%  vel: velocities on the vertices of the mesh
%  outdir: output directory (optional)
%
% Output:
%  strain: strain rates 
%
% Hua Wang @ Leeds, 17/11/2009
%============================================================
ntri=size(trim.tri,1);
incenter=zeros(ntri,2); %for gmt plot

for i=1:ntri
  %vertices id
  vtxid = trim.tri(i,:);
  %vertices coordinates for each triangle
  ll=[trim.x(vtxid),trim.y(vtxid)];

  %incenter of the triangular
  [incenter(i,:)]=tri2incenter(ll);

  %convert from lon/lat to x/y in km
  dxy=ll2utm(ll,ll(1,:));   %coordinates in km

  %calculate velocity gradient tensor
  %velgrd: [dvx/dx, dvx/dy, dvy/dx, dvy/dy]
  dvel = vel(vtxid',1:2)-repmat(vel(vtxid(1),1:2),3,1);
  dvel = reshape(dvel(2:3,:)',4,1);
  A=[dxy(2,:) 0 0; 0 0 dxy(2,:); dxy(3,:) 0 0; 0 0 dxy(3,:)];
  velgrad = A\dvel;

  %calculate strain rate and rotation rate tensor
  %deps: exx, eyy, exy=eyx, omega
  A=[1 0 0 0; 0 0 0 1; 0 0.5 0.5 0; 0 0.5 -0.5 0];
  deps(i,:)=A*velgrad;
end
deps=deps*1e-6;

%calculate principal strain rates
[peps]=prinstrain(deps);
strain=[deps,peps];

if nargin>2
  outstrain=[incenter,strain];
  save(strcat(outdir,'strain_tri.dat'),'outstrain','-ASCII');
end
