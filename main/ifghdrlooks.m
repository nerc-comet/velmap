function [ifghdr]=ifghdrlooks(ifghdr,lksx,lksy)
%===========================================
% update ifghdr after multi-look processing
%
% INPUT:
%   ifghdr: orignal ifghdr
%   lksx/y: looks number in x/y
% OUTPUT:
%   ifghdr: ifghdr after multi-look processing
%
% Hua Wang @ Uni Leeds, 04/11/2009
%===========================================

ifghdr.xfirst=ifghdr.xfirst+(lksx-1)*ifghdr.xstep/2;
ifghdr.yfirst=ifghdr.yfirst+(lksy-1)*ifghdr.ystep/2;
ifghdr.xstep = ifghdr.xstep*lksx;
ifghdr.ystep = ifghdr.ystep*lksy;
ifghdr.xpsize = ifghdr.xpsize*lksx;
ifghdr.ypsize = ifghdr.ypsize*lksy;
ifghdr.width = floor(ifghdr.width/lksx);
ifghdr.length = floor(ifghdr.length/lksy);
