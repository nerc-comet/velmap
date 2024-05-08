function [orbatm]=designorbatm(insar)
%==================================================================== 
%function [orbmat]=designorbatm(insar)
%                                                                     
% Make design matrix for orbit correction
%                                                                     
% INPUT:                                                              
%   insar:  insar data
%                                                                     
% OUTPUT:                                                             
%   orbatm: design matrix for orbital & atm errors                     
%                                                                     
% Hua Wang @ Uni Leeds, 23/09/2009
%
% 27/05/2011 HW: set orb/atm parameters one by one
%==================================================================== 

ninsar=length(insar);
orbatm=[];
for is=1:ninsar

  %-------------------------------
  % orbital errors design matrix
  %-------------------------------
  %make grid (x,y coordinates)
  [rows,cols]=size(insar(is).stackmap);
  n=rows*cols;
  if insar(is).proc.orbdegree==0
    B=ones(n,1);
  else
    [xx,yy]=meshgrid(1:cols,1:rows);
    xxv=reshape(xx',n,1);
    yyv=reshape(yy',n,1);
    %reference point
    refx=floor(cols/2);
    refy=floor(rows/2);
    xxv=xxv-refx;
    yyv=yyv-refy;
    xxv = xxv*insar(is).ifghdr.xpsize/100;  %test, using 100km as unit
    yyv = yyv*insar(is).ifghdr.ypsize/100;  %test, using 100km as unit
    clear('xx','yy');
 
    if insar(is).proc.orbdegree==1
      B=[xxv yyv ones(n,1)];
    else                    %degree == 2
      xxv2 = xxv.*xxv;
      yyv2 = yyv.*yyv;
      xyv  = xxv.*yyv;
      B=[xxv2 yyv2 xyv xxv yyv ones(n,1)];
    end
  end

  %-------------------------------
  % atm errors design matrix
  %-------------------------------
%  disp(insar(is).proc.atmdegree)
  if insar(is).proc.atmdegree==1
    demv=reshape(insar(is).dem',n,1);
    %disp(size(insar(is).dem))
    %disp(n)
    demv=demv/1000;  %using 1km as unit

    %remove nans
%    disp(size(B))
%    disp(size(demv))
    B(isnan(demv),:)=[];
    demv(isnan(demv),:)=[];
%    disp(size(B))
%    disp(size(demv))

    %remove reference height
    demv=demv-mean(demv);       

    iorbatm=[B demv];
  else
    ratev=reshape(insar(is).stackmap',n,1);
    B(isnan(ratev),:)=[];
    clear ratev;
    iorbatm=B;
  end

  %disp(size(iorbatm))
  orbatm = sparse(blkdiag(orbatm,double(iorbatm)));
%  disp(size(orbatm))
%  orbatm = sparse(orbatm);
%  disp(size(orbatm))
end
