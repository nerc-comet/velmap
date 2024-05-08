function []=writemat(filename,data,prec,interleave)
%==================================================
%function []=writemat(filename,data,prec,interleave)
%
% Write a matrix into a binary file    
%                                                  
% INPUT:                                           
%   filename: output file name                     
%   data: matrix to be written                   
%   prec: precision of input file, (optional, default 1: real4; 2: int16; 3: int8; others: double) 
%   interleave: the format in which the data is stored (default: 'real'; 'complex'; 'rmg')
% OUTPUT:                                          
%   NO                                             
%                                                  
% Hua Wang @ Uni Leeds, 09/07/2008                     
%                                                  
% 16/07/2011 HW: support more precisions and interleaves
% 27/08/2009 HW: support int16 format
%==================================================

fid = fopen(filename,'w','l');
if fid<0
  error([filename ' can not be opened!']);
end

if nargin<3
  prec=1;
end

[rows,cols,nbands]=size(data);
%sort out bands for the interleaves
if nargin<4 | nbands==1
  interleave='real';
end
interleave=lower(interleave);
if strcmp(interleave,'real')
  data_sort = data;
elseif strcmp(interleave,'rmg')
  data_sort=zeros(2*rows,cols);
  amp=[1:rows]*2-1;
  phs=[1:rows]*2;
  data_sort(amp,:)=data(:,:,1);
  data_sort(phs,:)=data(:,:,2);
else
  data_sort=zeros(rows,2*cols);
  amp=[1:cols]*2-1;
  phs=[1:cols]*2;
  data_sort(:,amp)=data(:,:,1);
  data_sort(:,phs)=data(:,:,2);
end

%output the data
if prec == 1
  sprec='real*4';
elseif prec == 2
  sprec='int16';
elseif prec==3
  sprec='int8';
else
  sprec='double';
end
fwrite(fid,data_sort',sprec);

fclose(fid);
