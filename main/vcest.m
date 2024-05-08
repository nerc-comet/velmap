function [ksqr,var0,var0_est,pb,pl,nbb,w]= vcest(res,x,x_vcm,var0_est,pb,pl,nbb,w)
%=================================================
%function [ksqr,var0,var0_est,pb,pl,nbb,w]= vcest(res,x,x_vcm,var0_est,pb,pl,nbb,w)
% Function on variance component estimation
% Inputs:
%   res: residuals
%   x: estimated unknown parameters
%   x_vcm: vcm of the unknown parameters
%   var0_est: inital estimated variance of unit weight
%   pb: inital weighted design matrix
%   pl: inital weighted obs
%   nbb: normal functions
%   w: constant
%
% Outputs:
%   ksqr: kai-square errors
%   var0: incremental variance of unit weight
%   var0_est: estimated variance of unit weight
%   pb: updated pb
%   pl: update pl
%   nbb: updated nbb
%   w: updated w
%
% Hua Wang, 30/09/2019, adopted from solve_vmp
%=================================================
ksqr=(pb*x-pl)'*res;
nobs=length(res);
r=nobs-trace(x_vcm*nbb);
var0=ksqr/r ;
var0_est=var0_est*var0; %resultant var0

%update pb, pl, nbb, w
pb=pb/var0;
pl=pl/var0;
nbb=nbb/var0;
w=w/var0;
