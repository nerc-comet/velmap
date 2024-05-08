
function [insar] = loadinsar(insarpar,procmode,tssmpar,gps)
%=============================================
%function [insar] = loadinsar(insarpar,procmode,tssmpar,gps)
%
% Load InSAR rate map and time series
%
% Input:
%   insarpar: insar parameters
%   procmode: processing mode (optional 1: static; 2: time series)
%   tssmpar: temporal smoothing parameters including t0/dt/mingps (optional)
%     t0 - intitial time of the time series
%     dt - time intervals of the time series
%     mingps - minimum number of GPS sites (default 5)
%   gps:  gps data (optional)
%
% Output:
%   insar: insar data structure, including stackmap etc.
%   insar=struct('proc',{},'ifghdr',{},'tsepoch',{},...
%               'vcm',{},'los',{},'azi',{},'dem',{});
%
% Hua Wang @ Uni Leeds, 21/09/2009
%
% 11/03/2015 HW: path structure for new pirate
% 02/09/2011 HW: load ts insar, and tidy by gps
% 27/05/2011 HW: set proc for each insar dataset
% 18/05/2011 HW: read rsc file instead of ifghdr.mat if available
%                using pixel size to determine looks number
% 07/10/2010 HW: save proper incidence file to 'out'
%=============================================
if nargin<2
  procmode=1;
end

if procmode==2 && nargin<3
  error('less than 3 input arguments for time series loading');
end

for i=1:insarpar.ninsarfile

  fprintf('preparing for the %d/%-d insar data ...\n',i,insarpar.ninsarfile);

  outdir=insarpar.dir(i);
  obsdir=char(strcat(outdir,'../obs/'));

  %read proc file for the insar rate processing, 27/05/2011 HW
  procfile=char(strcat(outdir,'insar.proc'));
  insar(i).proc=getinsarproc(procfile);
 
  %read header file
  %identify the version of pirate using the location of ifg.rsc
  ifgrscfile=char(strcat(outdir,'/ratemap/ifg.rsc'));  %pi-rate v3.2 and later on
  ifghdrfile=char(strcat(obsdir,'ifghdr.mat')); %old version of pi-rate 2.0 or earlier
  if exist(ifgrscfile,'file')
    ifghdr=rsc2hdr(ifgrscfile);
    ratemapdir=char(strcat(outdir,'/ratemap/'));
    auxdir=char(strcat(outdir,'/aux/'));
  else
    ifgrscfile=char(strcat(outdir,'ifg.rsc'));  %pi-rate before v3.2
    if exist(ifgrscfile,'file')
      ifghdr=rsc2hdr(ifgrscfile);
    elseif exist(ifghdrfile,'file') %pi-rate before v2.2 or so
      load(ifghdrfile);
    else
      error('ifg header file does not exist, please check ifghdr.mat or ifg.rsc');
    end
    ratemapdir=outdir;
    auxdir=outdir;
  end

  %determine looks by pixel size if available
  if (insarpar.xpsize~=0) & (insarpar.ypsize~=0)
    lksx=round(insarpar.xpsize/abs(ifghdr.xstep));
    lksy=round(insarpar.ypsize/abs(ifghdr.ystep));
  else
    lksx=insarpar.lksx;
    lksy=insarpar.lksy;
  end

  %load ifglist.mat to get vcm parameters
  ifglistfile=char(strcat(auxdir,'ifglist.mat'));
  if exist(ifglistfile)
    load(ifglistfile);
  end

  %-----------------------
  %  ratemap
  %-----------------------
  fprintf('loading rate map ...\n');
  stackmapname=char(strcat(ratemapdir,'stackmap.dat'));
  stackmap=readmat(stackmapname,ifghdr.length,ifghdr.width,1);
  insar(i).stackmap=looks(stackmap,lksx,lksy);
  clear('stackmap');

  %-----------------------
  %  incidence
  %-----------------------
  %prepare incidence file, modified on 07/12/2010
  %crop the incidence from a single file, 18/05/2011
  fprintf('loading incidence angle ... \n');
  incfile=char(strcat(ratemapdir,'incidence.unw'));
  %using incidence.unw
  if exist(incfile,'file')
    fileinfo=dir(incfile);
    filesize=fileinfo.bytes/8;
    if filesize==ifghdr.length*ifghdr.width
      inc=readmat(incfile,ifghdr.length,ifghdr.width,1,'rmg');
    else
      error('incidence filesize error (incidence.unw)');
    end
  %using a given file
  elseif exist(insar(i).proc.incfile,'file')
    inc=imagecrop(insar(i).proc.incfile,ifghdr,1,'rmg',0);
    multibandwrite(single(inc),incfile,'bil');
  %search all incidence files
  else
    inc=prepinc(obsdir,ifghdr,ifglist);
    multibandwrite(single(inc),incfile,'bil');
  end
  insar(i).los=looks(inc(:,:,1),lksx,lksy);
  %insar(i).azi=looks(inc(:,:,2),lksx,lksy);
  insar(i).azi=looks(inc(:,:,2) + 90,lksx,lksy); %JRW/HW hack
  insar(i).stackmap(isnan(insar(i).los))=nan;
  clear('inc');

  %-----------------------
  %  dem
  %-----------------------
  %read dem data
  %needn't dem file if atmdegree==0
  if insar(i).proc.atmdegree~=0
    fprintf('loading dem data ... \n');
    demname=char(strcat(ratemapdir,'dem.dat'));
    if exist(demname,'file')
      fileinfo=dir(demname);
      filesize=fileinfo.bytes/4;
      if filesize==ifghdr.length*ifghdr.width
        dem=readmat(demname,ifghdr.length,ifghdr.width,1);
      else
        error('dem filesize error (dem.dat)');
      end
    else
      dem=imagecrop(insar(i).proc.demfile,ifghdr,2);
      writemat(demname,dem);
    end
    insar(i).dem=looks(dem,lksx,lksy);
    clear('dem');
    insar(i).stackmap(isnan(insar(i).dem))=nan;
  end

  %-----------------------
  % time series
  %-----------------------
  if procmode==2

    %-----------------------
    % tsincr
    %-----------------------
    %get epochlist
    fprintf('loading time series ...\n');
    epochlistfile=char(strcat(auxdir,'epochlist.mat'));
    if exist(epochlistfile)
      load(epochlistfile);
    else
      error('must provide epochlist for time series mode');
    end
 
    %read time series
    nincr=size(epochlist.date,1)-1;
    tsdir=strcat(outdir,'timeseries/');
    tsincr=zeros(floor(ifghdr.length/lksy),floor(ifghdr.width/lksx),nincr);
    tserr=zeros(floor(ifghdr.length/lksy),floor(ifghdr.width/lksx),nincr);
    for its=1:nincr
      %read time series - tsincr
      itsfile=char(strcat(tsdir,'tsincr_',num2str(epochlist.date(its+1)),'.dat'));
      itsincr=readmat(itsfile,ifghdr.length,ifghdr.width,1);
      tsincr(:,:,its)=looks(itsincr,lksx,lksy);
      %read time series - tserror
      if insar(i).proc.errormap==0
        itsfile=char(strcat(tsdir,'tserr_',num2str(epochlist.date(its+1)),'.dat'));
        itserr=readmat(itsfile,ifghdr.length,ifghdr.width,1);
        tserr(:,:,its)=looks(itserr,lksx,lksy);
      end
    end

    %-----------------------
    %merge time series
    %-----------------------
    span=datenum(num2str(epochlist.date(2:nincr+1)),'yyyymmdd');
    breaks=tssmpar.t0-floor(tssmpar.dt/2):tssmpar.dt:span(nincr)+tssmpar.dt;
    [nhis,bin] = histc(span,breaks); %bin has the same size with intv
    pieces=length(breaks)-1;
    nhis(pieces)=nhis(pieces)+nhis(pieces+1);
    nhis(pieces+1)=[];
    breaks(pieces+1)=[];
    breaks=breaks+floor(tssmpar.dt/2);
    
    %set nhis=0 if no gps site exists
    if nargin>3
      gpsepoch=datenum(num2str(unique([gps.epoch]')),'yyyymmdd');
      nhis(ismember(breaks,gpsepoch)==0)=0;
    end

    %merge insar
    [yr,mm,dd]=datevec(breaks(nhis>0));
    insar(i).epoch=yr*10000+mm*100+dd;
    nep = length(insar(i).epoch);
    [nrows,ncols,~]=size(tserr);
    insar(i).tsvel=zeros(nrows,ncols,nep);
    mgerr=zeros(nrows,ncols,nep);
    intv=diff(epochlist.span);
    k=1;
    for j=1:pieces
      if nhis(j)>0
        index=find(bin==j);
        insar(i).tsvel(:,:,k)=sum(tsincr(:,:,index),3)/sum(intv(index));   %not accurate ??
        if insar(i).proc.errormap==0
          mgerr(:,:,k)=sqrt(sum(tserr(:,:,index).^2,3))/sum(intv(index))/10;   %not accurate ??
        end
        k=k+1;
      end
    end
    clear('tserr','tsincr');

    %---------------------------------
    %mask tsvel according to stackmap
    %---------------------------------
    %mask los/azi/dem by tsvel
    nincr=size(mgerr,3);
    tsmask=sum(insar(i).tsvel,3);
    insar(i).stackmap(isnan(tsmask))=nan;
    insar(i).tsvel(isnan(repmat(insar(i).stackmap,[1,1,nincr])))=nan;
    mgerr(isnan(insar(i).tsvel))=nan;

    %--------------------
    % time series error
    %--------------------
    %vcm, don't consider covariance between epochs ????
    fprintf('making vcm for the time series ...\n');
    insar(i).tsvcm=[];
    for its=1:nincr
      if insar(i).proc.errormap==0
        %%make vcm for each insar incremental
        %[maxvar,alpha]=make_vcm(insar(i).tsvel(:,:,its),lksx*ifghdr.xpsize,lksy*ifghdr.ypsize);
        %vcm_line=make_vcms(alpha,ifghdr.length,ifghdr.width,ifghdr.xpsize,ifghdr.ypsize,lksx,lksy,2);
        vcm_line=make_vcms(ifglist.alpha,ifghdr.length,ifghdr.width,ifghdr.xpsize,ifghdr.ypsize,lksx,lksy,2);
        vcm_its=make_vcms_ratemap(vcm_line,insar(i).tsvel(:,:,its),mgerr(:,:,its));
        insar(i).tsvcm=blkdiag(insar(i).tsvcm,vcm_its);
        clear('vcm_line','vcm_its');
      else                   %don't consider covariance
        tserr=insar(i).tsvel(:,:,its);
        tserr(~isnan(tserr))=insar(i).proc.errormap.^2;
        verr=reshape(tserr',numel(tserr),1);
        verr(isnan(verr))=[];
        insar(i).tsvcm=blkdiag(insar(i).tsvcm,sparse(double(diag(verr))));
        clear('tserr','verr');
      end
    end
  end

  %-----------------------
  %  errormap
  %-----------------------
  if insar(i).proc.errormap==0
    fprintf('loading error map ...\n');
    errormapname=char(strcat(ratemapdir,'errormap.dat'));
    errormap=readmat(errormapname,ifghdr.length,ifghdr.width,1);
    errormap=looks(errormap,lksx,lksy);
    errormap(isnan(insar(i).stackmap))=nan;
  else
    errormap=insar(i).stackmap;
    errormap(~isnan(errormap))=insar(i).proc.errormap;
  end
  %make vcm for each insar stackmap
  fprintf('making vcm for stackmap ... \n');  
 
 if exist(ifglistfile) %consider covariance
   fprintf('making vcm for stackmap ... considering covariance ... \n');
    %vcm_line=make_vcms(ifglist.alpha,ifghdr.length,ifghdr.width,ifghdr.xpsize,ifghdr.ypsize,lksx,lksy,2);
    %use mean alpha estimated using cvdcalc for off-fault frames - JW - Jul 2019
    %add last "3" n-times e-folding distance to speed up VCM create - JW - Aug 2019
    vcm_line=make_vcms(0.2,ifghdr.length,ifghdr.width,ifghdr.xpsize,ifghdr.ypsize,lksx,lksy,2);
    insar(i).vcm = make_vcms_ratemap(vcm_line,insar(i).stackmap,errormap);
 else %don't consider covariance
    fprintf('making vcm for stackmap ... not considering covariance ... \n');
    errormap=errormap.^2;
    verr=reshape(errormap',numel(errormap),1);
    verr(isnan(verr))=[];
    insar(i).vcm = sparse(double(diag(verr)));
    clear errormap,verr;
 end
 %added for parallel processing, HW
 insar(i).nobs=size(insar(i).vcm,1);

  %just for test, HW
  %insar(i).tsvcm=[];
  %for its=1:length(insar(i).epoch)
  %  insar(i).tsvcm=blkdiag(insar(i).tsvcm,insar(i).vcm);
  %end
    
  %update ifghdr
  insar(i).ifghdr=ifghdrlooks(ifghdr,lksx,lksy);
  insar(i).ifghdr
  %-----------------------
  %mask los/azi/dem by stackmap
  %-----------------------
  insar(i).los(isnan(insar(i).stackmap))=nan;
  insar(i).azi(isnan(insar(i).stackmap))=nan;
  if insar(i).proc.atmdegree~=0
    insar(i).dem(isnan(insar(i).stackmap))=nan;
  end
  end
end
