function [] = hdr2rsc(ifghdr,rscfile)
%==============================================
%function [] = hdr2rsc(ifghdr,rscfile)
%                                                                   
% Make ROI_PAC rsc file from a hdr matrix
%
% INPUT:
%   ifghdr: header parameters
%   rscfile: roi_pac rsc file
%
% OUTPUT:
%   None
%
% 21/11/2015 HW: output wavelength if existing
% Hua Wang, 12/04/2011
%==============================================

%only compatible for roi_pac format
%X_FIRST: topleft east
%Y_FIRST: topleft north

if nargin<2
  error('less than two input arguments');
end

fid=fopen(rscfile,'w');
if fid<0
  error(['can not open ' rscfile]);
end

fprintf(fid,'WIDTH             %-d\n',ifghdr.width);
fprintf(fid,'FILE_LENGTH       %-d\n',ifghdr.length);
fprintf(fid,'X_FIRST           %-16.9f\n',ifghdr.xfirst);
fprintf(fid,'X_STEP            %-16.9f\n',ifghdr.xstep);
fprintf(fid,'Y_FIRST           %-16.9f\n',ifghdr.yfirst);
fprintf(fid,'Y_STEP            %-16.9f\n',ifghdr.ystep);
if isfield(ifghdr,'wvl')
  fprintf(fid,'WAVELENGTH        %-16.9f\n',ifghdr.wvl);
end

fclose(fid);
