function[ifg_lowres]=looks(ifg,lksx,lksy,thrfrac)
%====================================================================
%function[ifg_lowres]=looks(ifg,lksx,lksy,thrfrac)
%                                                                    
% Multi-look for an input interferogram                     
%                                                                    
% INPUT:                                                             
%   ifg:  input interferogram                                        
%   lksx: looks number for x-component                               
%   lksy: looks number for y-component                               
%   thrfrac: fraction threshold (0-1) (default 0.5)
% OUTPUT:                                                            
%   ifg_lowres: multi-looked interferogram                           
%                                                                    
% Hua Wang @ Uni Leeds, 04/03/2008                                       
%                                                                    
% NOTE: threshold of valid pixels is set as half of the totle number 
%       in the window                                                
%
% 11/03/2015 HW: add threshold for valid pixels
% 27/04/2011 HW: set ifg_lowres as 'single' to save space
% 11/01/2011 HW: fix a bug while zeros are not NaNs
%====================================================================
if nargin<4
  thrfrac=0.3;
end
if thrfrac<0 || thrfrac>1
  error('thrfrac is wrong in looks: %f',thrfrac);
end
if (lksx==1 && lksy==1)
  ifg_lowres=ifg;
else
  [rows,cols]=size(ifg);
  rows_lowres=floor(rows/lksy);
  cols_lowres=floor(cols/lksx);
  ifg_lowres=NaN(rows_lowres,cols_lowres,'single');

  %threshold of valid pixels is set as half of the totle number in the window
  psize=lksx*lksy;
  thr = floor(psize*thrfrac);

  for r=1:rows_lowres
    for c=1:cols_lowres
      patch=ifg((r-1)*lksy+1:r*lksy,(c-1)*lksx+1:c*lksx);
      pv=reshape(patch,psize,1);
      n=sum(~isnan(pv));
      if n>=thr
        ifg_lowres(r,c)=nanmean(pv);
      end
    end
  end
end
