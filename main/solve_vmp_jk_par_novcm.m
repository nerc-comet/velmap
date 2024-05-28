function [fitmodel]=solve_vmp_jk_par_novcm(trim,smf,gps,insar,outdir,vce,njk,ncore)
%============================================================
%function [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smf,gps,insar,outdir,vce)
% 
% Design interpolation matrix for insar data
% 
% Input:
%  trim:      triangular mesh structure (tri,node)
%  smf:       smoothing factor
%  gps:       gps data
%  insar:     insar data structure
%  outdir:    output directory for the log file
%  vce:       variance component estimation (0: no vce; default 1)
%  njk:       number of jackknife samples
%
% Output:
%  fitmodel: fitted model, including velocity/orb/atm parameters
%  vcmmodel: covariance matrix for the model
%  wrss:   weighted residual squares sum
%  rough:  roughness
%
% Hua Wang @ Leeds, 23/09/2009
%
% 29/09/2015 HW: output resultant variance ration of InSAR vs GPS
% 12/03/2015 HW: variance component estimation is added
%============================================================

if nargin<6
  vce=1;
end

  tol=1.0e-8;
  maxit=500;
  fitmodel_local = [];
  ninsar=length(insar);
%disp(ninsar)
%pause
%delvce=1;
%iter=1;
%fid=fopen(strcat(outdir,'log'),'a');
%fprintf(fid,'---------smoothing factor: %f------------\n',smf);


%inverse parameters
invenu=getinvenu(gps,insar);

%vertex number
nvtx=length(trim.x);

%number of unknown parameters
npar=nvtx*sum(invenu);

%B=[];
%d=[];

ngf=length(gps);
%disp(ngf)

%%
disp('making design matrix for smoothing operator ...')
tic
[smmat]=designsmooth(trim,1,invenu);
%whos
smmat = sparse(smmat);
%whos
%pause
%size(smmat)
%disp(smmat(1:50,1:50))
%pause
smmat=smmat*smf;
%number of spatial smoothing observations for velocities
nsm=size(smmat,1); 
%ninsar=length(insar);
if ninsar>0
    orbatm=designorbatm(insar);
    npar_orbatm=size(orbatm,2);
    smmat=[smmat sparse(nsm,npar_orbatm)];
end
smmat = sparse(smmat);
%smoothing observation
nbb_sm=smmat'*smmat;
nbb_sm = sparse(nbb_sm);
%nsm = size(smmat,1);

  %B=[B;smmat];
  %B = sparse(B);
  %clear smmat
  %d=[d;zeros(nsm,1)];
  %whos
  %pause

%disp(ninsar)


%%
%--------------------------------------------
% weighted design matrix and obs of InSAR
%--------------------------------------------
if ninsar>0
  %design matrix for InSAR
  disp('making design matrix for insar data ...')
  parpool(8);
  insarmat=designinterp(trim,insar,0,invenu);
%  size(insarmat)

  %orb&atm error matrix
  orbatm=designorbatm(insar);
% disp(size(insarmat))
% disp(size(orbatm))

  %final design matrix for insar (sparse matrix)
  insarmat=[insarmat orbatm];
  %insarmat = sparse(insarmat);
%  size(insarmat)

  %insar observations and normal functions
  npar_orbatm=size(orbatm,2);
  npar=npar+npar_orbatm;
  nobsinsar=size(insarmat,1);
  disp('making observation vector and normal functions - a cup of coffee ...')
  i1=cumsum([insar.nobs]);
  i0=[1,i1(1:ninsar-1)+1];
  %myPool = parpool(4);
  parfor i=1:ninsar
%for i=1:ninsar
    fprintf('parallel processing %d/%d insar files ...\n',i,ninsar)
    vstackmap=double(reshape(insar(i).stackmap',[],1));
    vstackmap(isnan(vstackmap))=[];
    %tmp(i).vinsar=vstackmap; % tmp for parallel computing
 
    %weighted design matrix and obs
    tmp(i).pb=choldiv(insar(i).vcm,insarmat(i0(i):i1(i),:));
    tmp(i).pl=choldiv(insar(i).vcm,vstackmap);
  end
  %delete(gcp('nocreate'))
  %vinsar=double(vertcat(tmp.vinsar));
  pb_insar=vertcat(tmp.pb);
  pl_insar=vertcat(tmp.pl);
  clear('tmp','vstackmap');
  nbb_insar=insarmat'*pb_insar;
  w_insar=insarmat'*pl_insar;
  %B=[B;insarmat];
  %d=[d;vinsar];
  clear insarmat orbatm insar pb_insar pl_insar
else
  nobsinsar=0;
  nbb_insar = [];
  w_insar = [];
end
delete(gcp('nocreate'))
%return
%w=sparse(npar,1);
%B = sparse(B);
%whos
%pause

%%
%--------------------------------------------
% weighted design matrix and obs of GPS
%--------------------------------------------
disp('making design matrix and obs of gps data ...')
  i1=cumsum([gps.nsite].*[gps.ndim]);
  i0=[1,i1(1:ngf-1)+1];
  nobsgps=i1(ngf);
  vgps=zeros(nobsgps,1);
  gpsmat=[];
  gpsvcm=[];
  for i=1:ngf
    %obs
    vgps(i0(i):i1(i))=reshape(vertcat(gps(i).site.vel),[],1);
    %disp(i0(i))
    %disp(i1(i))
 
    %design matrix
    igpsmat=designgps(trim,gps(i).site,gps(i).invenu);
    %append zeros for 2D gps files
    if (invenu(3)==1 && gps(i).invenu(3)==0)
      igpsmat=[igpsmat,sparse(size(igpsmat,1),nvtx)];
    end
    if ninsar>0
      igpsmat=[igpsmat sparse(size(igpsmat,1),npar_orbatm)];
    end
    gpsmat=[gpsmat;igpsmat];
 
    %vcm
    igpsvcm=sparse(size(igpsmat,1),size(igpsmat,1));
    for is=1:gps(i).nsite
      index=[is:gps(i).nsite:is+(gps(i).ndim-1)*gps(i).nsite];
      igpsvcm(index',index)=gps(i).site(is).vcm;
    end
    gpsvcm=blkdiag(gpsvcm,igpsvcm);
    clear('index','igpsmat','igpsvcm');
  end
  %clear gps
  %Nbb and W
  gpsvcm = sparse(gpsvcm);
  gpsmat = sparse(gpsmat);

  %pb_gps=choldiv(gpsvcm,gpsmat);
  %pl_gps=choldiv(gpsvcm,vgps);
%  disp(vgps)
%  disp(vgps(1:100))
%  pause
%  gpssize = size(gpsvcm,1);
  bootindices = cell(njk,1);
  whos
  %pause
  %parpool(4);
      %rng('shuffle');
  for jkindex = 1:njk
      rng('shuffle');
      bootindices{jkindex} = [];
      %bootindices = [];
      for i = 1:ngf
          %rng('shuffle');
          %disp(gps(i).nsite)
          bootlength = round(rand(1,1)*(gps(i).nsite - 1)) + 1;
          while bootlength<10
          bootlength = round(rand(1,1)*(gps(i).nsite - 1)) + 1;
          end
          stationindices = sort(randperm(gps(i).nsite,bootlength))';
          %disp(stationindices)
          if gps(i).ndim==2
              bootindices_local = sort([stationindices*2 - 1; stationindices*2]);
          elseif gps(i).ndim==3
              bootindices_local = sort([stationindices*3 - 2; stationindices*3 - 1; stationindices*3]);
          end
          if i>1
              %disp(gps(1).nsite*gps(1).ndim)
              bootindices_local = bootindices_local + gps(1).nsite*gps(1).ndim;
          end
          bootindices{jkindex} = [bootindices{jkindex}; bootindices_local];
      end
  end

  %disp(size(gpsmat))
  %disp(size(gpsvcm))
  %disp(size(vgps))

  parpool(ncore);
  parfor jkindex = 1:njk
    fprintf('jkindex %d\n',jkindex)

    %disp(bootindices)
      %pause
      
      gpsmat_local = gpsmat(bootindices{jkindex},:);
      gpsvcm_local = gpsvcm(bootindices{jkindex},bootindices{jkindex});
      gpsvcm_sqrtdiag = sqrt(diag(gpsvcm_local));
      disp(max(gpsvcm_sqrtdiag))
      gpsvcm_sqrtdiag(gpsvcm_sqrtdiag>1e5) = 0;
      disp(max(gpsvcm_sqrtdiag))

      vgps_local = vgps(bootindices{jkindex});
      disp(vgps_local(1:10))
      vgps_local = vgps_local + randn(length(vgps_local),1).*gpsvcm_sqrtdiag;
      disp(vgps_local(1:10))

      pb_gps_local = choldiv(gpsvcm_local,gpsmat_local);
      pl_gps_local = choldiv(gpsvcm_local,vgps_local);
      nbb_local = sparse(gpsmat_local'*pb_gps_local);
      w_local = gpsmat_local'*pl_gps_local;
      gpsmat_local = [];
      gpsvcm_local = [];
      vgps_local = [];
      pb_gps_local = [];
      pl_gps_local = [];
      % bootlength = [];
      % bootindices = [];
      %clear gpsmat_local gpsvcm_local vgps_local pb_gps_local pl_gps_local bootlength bootindices

      nbb_local=nbb_sm+nbb_local;

      %fprintf('ninsar is %f\n',ninsar)
      if ninsar>0
          nbb_local=nbb_local+nbb_insar;
          w_local=w_local+w_insar;
      end
      nbb_local = sparse(nbb_local);

      %     disp(size(nbb))
      %     disp(size(w))
      % pause
      %--------------------------------------------
      % solve the system of equations
      %--------------------------------------------
      %whos
        fprintf('solving l1 and u1 in jk %d\n',jkindex)
        [l1,u1]=ilu(nbb_local,struct('type','ilutp','droptol',tol));
        fprintf('solving model in jk %d\n',jkindex)
     fitmodel_local(:,jkindex)=bicg(nbb_local,w_local,tol,maxit,l1,u1);
%      clear l1 u1 nbb_local w_local nbb_gps_local w_gps_local
  end
       clear l1 u1 nbb_local w_local nbb_gps_local w_gps_local
 
       %disp(size(fitmodel_local))
       disp(fitmodel_local(1:1000:end,:))
       disp(size(fitmodel_local))
       fitmodel_local = [unique(fitmodel_local','rows')]';
       disp(size(fitmodel_local))
%save('fitmodel_local.mat','fitmodel_local','-v7.3')
  save(strcat(outdir,'fitmodel_local.mat'),'fitmodel_local','-v7.3');

  sel = max(abs(fitmodel_local),[],1)<=45;
disp(sum(sel))
disp(size(fitmodel_local))
disp(size(fitmodel_local(:,sel)))
fitmodel = mean(fitmodel_local(:,sel),2);

  %fitmodel = median(fitmodel_local,2);
  %fitmodel = mean(fitmodel_local,2);
  clear fitmodel_local
  %delete(gcp('nocreate'))

  % nbb_gps=gpsmat'*pb_gps;
  % w_gps=gpsmat'*pl_gps;
  % %B=[B;gpsmat];
  % %d=[d;vgps];
  % 
  % nbb=nbb_sm+nbb_gps;
  %   clear nbb_sm
  %   %w=w+w_gps;
  % if ninsar>0
  %   nbb=nbb+nbb_insar;
  %   %w=w+w_insar;
  % end
  % nbb = sparse(nbb);

  %temptest = flip(sort(abs(fitmodel)));
  %disp(temptest(1:500))

  %vcm model
%  disp('generating vcm of model ...')
%  vcmmodel=choldiv(nbb,speye(size(nbb,1)));
%  vcmmodel = sparse(vcmmodel);

% disp('making design matrix for smoothing operator ...')
% tic
% [smmat]=designsmooth(trim,1,invenu);
% %size(smmat)
% %disp(smmat(1:50,1:50))
% %pause
% smmat=smmat*smf;
% %number of spatial smoothing observations for velocities
% nsm=size(smmat,1); 
% if ninsar>0
%   smmat=[smmat sparse(nsm,npar_orbatm)];
% end
% %smoothing observation
% vsm=zeros(nsm,1);
% 
   %B=[B;smmat];
   %B = sparse(B);
%   clear smmat
   %d=[d;zeros(nsm,1)];

%number of geodetic obs
%m=nobsinsar+nobsgps;
%total observation number
%nobs=m+nsm;


%    nbb=nbb_sm;
  %w=sparse(npar,1);
    % disp(size(nbb))
    % disp(size(w))
%   opts.POSDEF = true;
%   opts.SYM = true;
%   vcmmodel = linsolve(full(nbb),eye(size(nbb)),opts);

  %residuals
%   disp('generating weighted residuals ...')
%   res = B*fitmodel-d;
% 
%   %--------------------------------------------
%   % reweighting
%   %--------------------------------------------
%   %unit variance for GPS
%   disp('variance components estimating ...')
% 
%   var0_insar_est=1;
% var0_gps_est=1;
% 
%   res_gps=res(1+nobsinsar:nobsinsar+nobsgps);
%   [ksqr_gps,var0_gps,var0_gps_est,~,~,~,~]...
%    = vcest(res_gps,fitmodel,vcmmodel,var0_gps_est,pb_gps,pl_gps,nbb_gps,w_gps);
%   fprintf(fid,'current and cumulative GPS postpriori variance of unit weight: %5.2f, %5.2f\n',var0_gps,var0_gps_est);
% 
% %       disp(size(var0_gps))
% %       disp(size(vcmmodel))
% %  var0_gps = sparse(var0_gps);
% %  vcmmodel = sparse(vcmmodel);
% %         disp(size(var0_gps))
% %       disp(vcmmodel)
% 
%   %delvce=0;
%   ksqr_insar=0;
%   if ninsar>0
%     res_insar=res(1:nobsinsar);
%   end
%   if vce==0 || ninsar==0
%       fprintf('reweighting using estimated variance of unit weight');
% %      whos
%     %vcmmodel=var0_gps*vcmmodel; %reweight using estimated variance of unit weight
% %        whos
%   else
%     %unit variance for InSAR
%     if ninsar>0
%       [ksqr_insar,var0_insar,var0_insar_est,~,~,~,~]...
%        = vcest(res_insar,fitmodel,vcmmodel,var0_insar_est,pb_insar,pl_insar,nbb_insar,w_insar);
% %      delvce=delvce+abs(var0_gps/var0_insar-1);
%       fprintf(fid,'current and cumulative InSAR postpriori variance of unit weight: %5.2f, %5.2f\n',var0_insar,var0_insar_est);
%     end
% %    fprintf(fid,'variance component estimation - interation=%3d, delvce=%4.2f\n',iter,delvce);
% 
% %    iter=iter+1;
%   end
% %end
% 
% %--------------------------------------------
% % output RMS vs Roughness
% %--------------------------------------------
% %RMS GPS
% if ngf>0
%   rms_gps=sqrt(mean(res_gps.^2));
% %  i1=cumsum([gps.nsite].*[gps.ndim]);
%   i0=[1,i1(1:ngf-1)+1];
%   for i=1:ngf
%     ires=reshape(res_gps(i0(i):i1(i)),[],gps(i).ndim);
%     rms_gps_enu=sqrt(mean(ires.^2));
%     if gps(i).ndim==2
%       rms_gps_enu(3)=nan;
%     end
%     fprintf(fid,'RMS_GPS_F%d (E, N, U, Total): %f %f %f %f\n',i,rms_gps_enu,rms_gps);
%   end
% end
% 
% %RMS InSAR
% if ninsar>0
%   rms_insar=sqrt(mean(res_insar.^2));
%   fprintf(fid,'RMS_InSAR: %f\n',rms_insar);
% end
% 
% %RMS
% rms_all=sqrt(mean(res(1:m).^2));
% fprintf(fid,'RMS ALL: %f\n',rms_all);
% 
% %wrss & roughness
% if nargout>2
%   wrss=sqrt((ksqr_gps+ksqr_insar)/m);
%   rough=sqrt(res(m+1:nobs)'*res(m+1:nobs)/nsm)/smf;
%   fprintf(fid,'weighted rss: %f\n',wrss);
%   fprintf(fid,'roughness: %f\n',rough);
% end

%fclose(fid);
