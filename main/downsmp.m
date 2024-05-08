function [smp,ix,iy]=downsmp(img,lksx,lksy)
%=============================================
%function [smp,ix,iy]=downsmp(img,lksx,lksy)
%
% Down-sample a image
%
% Input: 
%  img:  input image with full resolution
%  lksx: looks number in x
%  lksy: looks number in y
%
% Output:
%  smp:  down-sampled image
%  ix/iy: donw-sampled pixel index
%
% Hua Wang @ Uni Leeds, 21/09/2009
%=============================================
[rows,cols]=size(img);
sx=floor(lksx/2);
sy=floor(lksx/2);
ix=(1:floor(cols/lksx))*lksx-sx;
iy=(1:floor(rows/lksy))*lksy-sy;
smp=img(iy,ix);
