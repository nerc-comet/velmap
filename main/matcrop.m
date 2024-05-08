function [out] = matcrop(inmat,inhdr,outhdr)
%===================================================================
%function [out] = matcrop(inmat,inhdr,outhdr)
%                                                                   
% Crop/Fill an image to the demension defined by a header file
%
% INPUT:
%   inmat:  input image (matrix)
%   inhdr:  input region
%   outhdr: target region
%
% OUTPUT:
%   out:  output data
%
% Hua Wang, 08/05/2013
%===================================================================

%import parameters in the resource file for each data
%only compatible for roi_pac format
%X_FIRST: topleft east
%Y_FIRST: topleft north
inhdr.xlast=inhdr.xfirst+(inhdr.width-1)*inhdr.xstep;
inhdr.ylast=inhdr.yfirst+(inhdr.length-1)*inhdr.ystep;

lksx = round(outhdr.xstep/inhdr.xstep);
lksy = round(outhdr.ystep/inhdr.ystep);
xlast=outhdr.xfirst+(outhdr.width-1)*outhdr.xstep;
ylast=outhdr.yfirst+(outhdr.length-1)*outhdr.ystep;

%get new dimensions for the output data
%note: using original pixel spacing, otherwise it may be wrong
%width/length in original resolution
orgwidth=outhdr.width*lksx;
orglength=outhdr.length*lksy;
dnx0 = floor((outhdr.xfirst-inhdr.xfirst)/inhdr.xstep);
if dnx0>0
  ix0 = dnx0+1;
  x0 = 1;
else
  ix0 = 1;
  x0 = -dnx0+1;
end

dny0 = floor((outhdr.yfirst-inhdr.yfirst)/inhdr.ystep);
if dny0>0
  iy0 = dny0+1;
  y0 = 1;
else
  iy0 = 1;
  y0 = -dny0+1;
end

dnx1 = floor((xlast-inhdr.xlast)/inhdr.xstep);
if dnx1>0
  ix1 = inhdr.width;
  x1 = orgwidth-dnx1;
else
  ix1 = inhdr.width+dnx1;
  x1 = orgwidth;
end

dny1=floor((ylast-inhdr.ylast)/inhdr.ystep);
if dny1>0
  iy1 = inhdr.length;
  y1 = orglength - dny1;
else
  iy1 = inhdr.length + dny1;
  y1 = orglength;
end

%it will ocassionally go wrong due to multi-look processing
%update nx,ny will make the program more stable
nx = min(x1-x0,ix1-ix0);
ny = min(y1-y0,iy1-iy0);
x1 = x0+nx;
y1 = y0+ny;
ix1 = ix0+nx;
iy1 = iy0+ny;

%crop in original resolution
nband=size(inmat,3);
croporg=NaN(orglength,orgwidth,nband,'single');
croporg(y0:y1,x0:x1,:)=inmat(iy0:iy1,ix0:ix1,:);
clear inmat;

%multilook
for i=1:nband
  out(:,:,i)=looks(croporg(:,:,i),lksx,lksy);
end
