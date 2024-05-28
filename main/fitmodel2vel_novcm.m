function [fitvtx] = fitmodel2vel_novcm(trim,fitmodel,invenu,outdir)
%========================================================
%function [fitvtx] = fitmodel2vel(trim,fitmodel,vcmmodel,invenu,outdir)
%
%  convert from fitted model to velocity field
%
% INPUT:
%  trim:     triangular mesh
%  fitmodel: fitted velocity field
%  vcmmodel: vcm of fitted velocity field
%  invenu:   inversion parameters
%  outdir:   output directory (optional)
% 
% OUTPUT:
%  fitvtx:  fitted velocity field, same format with gps
%
% Hua Wang @ Uni Leeds, 09/11/2009
%
% 27/05/2011 HW: deal with 3d velocities
%========================================================

%velocities on the vertics
nvtx=length(trim.x);
invs=sum(invenu);
vel=reshape(fitmodel(1:invs*nvtx),nvtx,invs);

% %error bar
% sigmavel=zeros(nvtx,invs*(invs+1)/2);
% for i=1:nvtx
%   %std
%   for j=1:invs
%     sigmavel(i,j)=sqrt(vcmmodel(i+(j-1)*nvtx,i+(j-1)*nvtx));      %stdx
%   end
%   %cov
%   l=invs+1;
%   for j=1:invs-1
%     for k=j+1:invs
%        sigmavel(i,l)=vcmmodel(i+(j-1)*nvtx,i+(k-1)*nvtx)/sigmavel(i,j)/sigmavel(i,k); %covariance coefficient
%        l=l+1;
%     end
%   end
% end

if nargout>0
  %3d velocities and uncertainties
  vel3d=zeros(nvtx,3);
  vel3d(:,invenu~=0)=vel;
  %sigmavel3d=zeros(nvtx,6);
  %sigpar = [invenu,invenu(1)*invenu(2),invenu(1)*invenu(3),invenu(2)*invenu(3)];
  %sigmavel3d(:,sigpar~=0)=sigmavel;

  %temporally use
  for i=1:nvtx
    fitvtx(i).lon=trim.x(i);
    fitvtx(i).lat=trim.y(i);
    fitvtx(i).vel=vel3d(i,:);
    %fitvtx(i).vel=vel(i,:);
    fitvtx(i).staid=num2str(i);
  end
end

if nargin>3
  %save velocity field data as the same format with gps input file
  %can use loadgps.m to read this file
  vtxid = [1:nvtx]';
  grdvel=[trim.x trim.y vel 0*vel vtxid];
  save(strcat(outdir,'velfit.dat'),'grdvel','-ASCII');

  %keep in matlab format for orb/atm etc.
  save(strcat(outdir,'fitmodel.mat'),'fitmodel');
  %save(strcat(outdir,'vcmmodel.mat'),'vcmmodel','-v7.3');
end
