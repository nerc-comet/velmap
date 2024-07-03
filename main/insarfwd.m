function [insarfit] = insarfwd(insar,trim,fitmodel,invenu,outdir,gps)
%========================================================
%function insarfwd(insar,trim,fitmodel,invenu,outdir,gps)
%
%  plot & output fitted insar data
%
% INPUT:
%  insar:     insar observations
%  trim:      triangular mesh
%  fitmodel:  fitted velocity field
%  invenu:    inversion parameters
%  outdir:    output directory
%  gps:       gps observations (optional)
%
% OUTPUT:
%  insarfit:  fitted insar data
%
% Hua Wang @ Uni Leeds, 09/11/2009
%
% 14/03/2014 HW: add comparison between InSAR and GPS
%========================================================

if nargin>5
  ngf=length(gps);
else
  ngf=0;
end

nvtx=length(trim.x);
invs=sum(invenu);
vel_v=fitmodel(1:invs*nvtx);

%insar interpolation matrix
insarmat=designinterp(trim,insar,0,invenu);
vinsarfit=insarmat*vel_v;

j=0;           %start position of observations
k=invs*nvtx;   %start position of orb/atm parameters
ninsar=length(insar);

for i=1:ninsar

  norbcoef=(insar(i).proc.orbdegree+1)*(insar(i).proc.orbdegree+2)/2;

  %find valid pixels
  [row,col]=size(insar(i).stackmap);
  vstackmap=reshape(insar(i).stackmap',row*col,1);
  idx=find(~isnan(vstackmap));

  %make folder
  ioutdir=char(strcat(outdir,'/insarfit',num2str(i,'%02d'),'/'));
  if ~exist(ioutdir,'dir')
    mkdir(ioutdir)
  end

  %ratemap
  nvalid=length(idx);  %number observations for the current insar stackmap
  vratemap=nan(row*col,1);
  vratemap(idx)=vinsarfit(j+1:j+nvalid);
  insarfit(i).ratemap=reshape(vratemap,col,row)';
  writemat(char(strcat(ioutdir,'ratemap.dat')),insarfit(i).ratemap);

  %orbmap
  orbatm=designorbatm(insar(i));
  orbmat=orbatm(:,1:norbcoef);
  orbpar=fitmodel(k+1:k+norbcoef);
  vorbmap=nan(row*col,1);
  vorbmap(idx)=orbmat*orbpar;
  insarfit(i).orbmap=reshape(vorbmap,col,row)';
  writemat(char(strcat(ioutdir,'orbmap.dat')),insarfit(i).orbmap);
  k=k+norbcoef;  %update position of orb/atm position

  %atmmap
  if insar(i).proc.atmdegree~=0
    atmpar=fitmodel(k+1);
    atmmat=orbatm(:,norbcoef+1);
    vatmmap=nan(row*col,1);
    vatmmap(idx)=atmmat*atmpar;
    insarfit(i).atmmap=reshape(vatmmap,col,row)';
    writemat(char(strcat(ioutdir,'atmmap.dat')),insarfit(i).atmmap);
    k=k+1;
  end

  %stackmap
  insarfit(i).stackmap=insarfit(i).ratemap+insarfit(i).orbmap;
  if insar(i).proc.atmdegree~=0
    insarfit(i).stackmap=insarfit(i).stackmap+insarfit(i).atmmap;
  end
  writemat(char(strcat(ioutdir,'stackmap.dat')),insarfit(i).stackmap);

  %residuals
  insarfit(i).resmap=insar(i).stackmap-insarfit(i).stackmap;
  writemat(char(strcat(ioutdir,'resmap.dat')),insarfit(i).resmap);

  %obs.
  writemat(char(strcat(ioutdir,'stackobs.dat')),insar(i).stackmap);

  %output header
  insarfit(i).ifghdr=insar(i).ifghdr;
  hdr2rsc(insar(i).ifghdr,char(strcat(ioutdir,'ifg.rsc')));

  %-----------------------
  %compare InSAR and GPS
  %-----------------------
  fid=fopen(char(strcat(ioutdir,'gpsvsinsar.dat')),'w');
  for igf=1:ngf
    igps=gps(igf).site;
    gpsdim=size(igps(1).vel,2);
    %remove gps points outside track
    x = floor(([igps.lon]-insar(i).ifghdr.xfirst)/insar(i).ifghdr.xstep) + 1;
    y = floor(([igps.lat]-insar(i).ifghdr.yfirst)/insar(i).ifghdr.ystep) + 1;
    tidy = find(x<1 | x>insar(i).ifghdr.width | y<1 | y>insar(i).ifghdr.length);
    gpsinside = igps;
    gpsinside(tidy)=[];
    x(tidy)=[];
    y(tidy)=[];
    %remove gps points without corresponding insar obs
    ngpsleft=length(gpsinside);
    for ig = 1:ngpsleft
      ipatxmin = max(1,x(ig)-1);
      ipatymin = max(1,y(ig)-1);
      ipatxmax = min(insar(i).ifghdr.width,x(ig)+1);
      ipatymax = min(insar(i).ifghdr.length,y(ig)+1);
      absratemap=insarfit(i).ratemap+insarfit(i).resmap;
      ipatsar = absratemap(ipatymin:ipatymax,ipatxmin:ipatxmax);
      if sum(sum(isnan(ipatsar)))<numel(ipatsar)
        ipatlos = insar(i).los(ipatymin:ipatymax,ipatxmin:ipatxmax);
        mlos=mean(mean(ipatlos(~isnan(ipatsar))));
        ipatazi = insar(i).azi(ipatymin:ipatymax,ipatxmin:ipatxmax);
        mazi=mean(mean(ipatazi(~isnan(ipatsar))));
        msar=nanmean(nanmean(ipatsar));
        %project gps to los
        uvec=unitvec(mlos,mazi);
        if gpsdim==2%remove Up component for 2d GPS
          uvec(3)=[];
        end
        gpslos=gpsinside(ig).vel*uvec';
        gpslos_sigma=sqrt(uvec*gpsinside(ig).vcm*uvec');
        fprintf(fid,'%f %f %f %f %f %s\n',msar,gpsinside(ig).lon,gpsinside(ig).lat,gpslos,gpslos_sigma,char(gpsinside(ig).staid));
      end
    end
  end
  fclose(fid);

  j=j+nvalid;
end
