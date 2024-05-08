function [fitmodel,vcmmodel,wrss,rough]=solve_vmp(trim,smf,gps,insar,outdir,invenu,vce)
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
%  invenu:    components to invert for (east, north, up)
%  vce:       variance component estimation (0: no vce; default 1)
%
% Output:
%  fitmodel: fitted model, including velocity/orb/atm parameters
%  vcmmodel: covariance matrix for the model
%  wrss:   weighted residual squares sum
%  rough:  roughness
%
% Hua Wang @ Leeds, 23/09/2009
%
% 10/08/2021 AW: added invenu as an imput
% 29/09/2015 HW: output resultant variance ration of InSAR vs GPS
% 12/03/2015 HW: variance component estimation is added
%============================================================

if nargin<7
    vce=1;
end

%inverse parameters
invenu=getinvenu(gps,insar);

%vertex number
nvtx=length(trim.x);

%number of unknown parameters
npar=nvtx*sum(invenu);

B=[];
d=[];
%--------------------------------------------
% weighted design matrix and obs of InSAR
%--------------------------------------------
ninsar=length(insar);
if ninsar>0
    
    % design matrix for InSAR
    disp('making design matrix for insar data ...')
    insarmat=designinterp(trim,insar,0,invenu);
    
    % orb & atm error matrix
    orbatm=designorbatm(insar);
    
    % final design matrix for insar (sparse matrix)
    insarmat=[insarmat orbatm];
    
    % insar observations and normal functions
    npar_orbatm = size(orbatm,2);
    npar = npar + npar_orbatm;
    nobsinsar = size(insarmat,1);
    disp('making observation vector and normal functions - a cup of coffee ...')
    i1 = cumsum([insar.nobs]);
    i0 = [1,i1(1:ninsar-1)+1];
    
    parfor i=1:ninsar
        % obs
        fprintf('parallel processing %d/%d insar files ...\n',i,ninsar)
        vstackmap = double(reshape(insar(i).stackmap',[],1)); % insar vels in column
        vstackmap(isnan(vstackmap)) = [];
        tmp(i).vinsar = vstackmap; % tmp for parallel computing
        
        % weighted design matrix and obs
        tmp(i).pb = choldiv(insar(i).vcm,insarmat(i0(i):i1(i),:));
        tmp(i).pl = choldiv(insar(i).vcm,vstackmap);
    end
    
    vinsar = double(vertcat(tmp.vinsar));
    pb_insar = vertcat(tmp.pb);
    pl_insar = vertcat(tmp.pl);
    clear('tmp','vstackmap');
    nbb_insar = insarmat'*pb_insar;
    w_insar = insarmat'*pl_insar;
    B = [B;insarmat];
    d = [d;vinsar]; % all the insar vels
else
    nobsinsar=0;
end

%--------------------------------------------
% weighted design matrix and obs of GPS
%--------------------------------------------
disp('making design matrix and obs of gps data ...')
ngf=length(gps);
if ngf>0
    i1=cumsum([gps.nsite].*[gps.ndim]);
    i0=[1,i1(1:ngf-1)+1];
    nobsgps=i1(ngf);
    vgps=zeros(nobsgps,1);
    gpsmat=[];
    gpsvcm=[];
    for i=1:ngf
        % obs
        vgps(i0(i):i1(i))=reshape(vertcat(gps(i).site.vel),[],1);           
        
        % design matrix
        igpsmat=designgps(trim,gps(i).site,gps(i).invenu);
        
        % append zeros for 2D gps files
        if invenu(3)==1 && gps(i).invenu(3)==0
            igpsmat=[igpsmat sparse(size(igpsmat,1),nvtx)];
%             igpsmat=[igpsmat; sparse(gps(i).nsite,size(igpsmat,2))];
        end
        
        % remove components that we don't wish to invert for
%         cind2 = cumsum(repelem(nvtx,3)); %cols
%         cind1 = [0 cind2(1:2)] + 1;
%         rind2 = cumsum(repelem(gps(i).nsite,3)); % rows
%         rind1 = [0 rind2(1:2)] + 1;
%         for jj = 1:3 % for each component
%             if invenu(jj) == 0
%                 igpsmat(rind1(jj):rind2(jj),:) = [];
%                 igpsmat(:,cind1(jj):cind2(jj)) = [];
%             end
%         end
        
        % append zeros for orbit and atmo
        if ninsar>0
            igpsmat=[igpsmat sparse(size(igpsmat,1),npar_orbatm)];
        end
        
        gpsmat=[gpsmat;igpsmat];
        
        % vcm
        igpsvcm=sparse(size(igpsmat,1),size(igpsmat,1));
        for is=1:gps(i).nsite
            index = is : gps(i).nsite : is+(gps(i).ndim-1)*gps(i).nsite;
            igpsvcm(index',index) = gps(i).site(is).vcm;
        end
        gpsvcm=blkdiag(gpsvcm,igpsvcm);
        clear('index','igpsmat','igpsvcm');
    end
    
    % Nbb and W
    pb_gps = choldiv(gpsvcm,gpsmat);
    pl_gps = choldiv(gpsvcm,vgps);
    nbb_gps = gpsmat'*pb_gps;
    w_gps = gpsmat'*pl_gps;
    B = [B;gpsmat];
    d = [d;vgps];
    clear('gpsvcm');
else
    nobsgps=0;
end

%-------------------------------------
%spatial smoothing matrix and obs
%-------------------------------------
disp('making design matrix for smoothing operator ...')
[smmat]=designsmooth(trim,1,invenu);
smmat=smmat*smf;

%number of spatial smoothing observations for velocities
nsm=size(smmat,1);
if ninsar>0
    smmat=[smmat sparse(nsm,npar_orbatm)];
end

%smoothing observation
vsm=zeros(nsm,1);
nbb_sm=smmat'*smmat;
B=[B;smmat];
d=[d;vsm];

%number of geodetic obs
m=nobsinsar+nobsgps;

%total observation number
nobs=m+nsm;

%-------------------------------------
% solve the system of equations
%-------------------------------------
var0_insar_est = 1;
var0_gps_est = 1;
delvce = 1;
iter = 1;
fid = fopen(strcat(outdir,'log'),'a');
fprintf(fid,'---------smoothing factor: %f------------\n',smf);
while (delvce>0.1 && iter<10)
    %--------------------------------------------
    % combine all normal functions
    %--------------------------------------------
    nbb=nbb_sm;
    w=sparse(npar,1);
    if ngf>0
        nbb = nbb+nbb_gps;
        w = w+w_gps;
    end
    if ninsar>0
        nbb = nbb+nbb_insar;
        w = w+w_insar;
    end
    
    %--------------------------------------------
    % solve the system of equations
    %--------------------------------------------
    disp('solving the system of equations ...')
    tol=1.0e-8;
    maxit=500;
    nbb = sparse(nbb);
    [l1,u1] = ilu(nbb,struct('type','ilutp','droptol',tol));
    fitmodel = bicg(nbb,w,tol,maxit,l1,u1); 
    clear l1 u1;
    
    % vcm model
    disp('generating vcm of model ...')
    vcmmodel = choldiv(nbb,speye(size(nbb,1)));
    %   opts.POSDEF = true;
    %   opts.SYM = true;
    %   vcmmodel = linsolve(full(nbb),eye(size(nbb)),opts);
    
    % residuals
    disp('generating weighted residuals ...')
    res = B*fitmodel-d;
    
    %--------------------------------------------
    % reweighting
    %--------------------------------------------
    % unit variance for GPS
    disp('estimating variance components ...')
    
    res_gps = res(1+nobsinsar:nobsinsar+nobsgps);
    [ksqr_gps,var0_gps,var0_gps_est,pb_gps,pl_gps,nbb_gps,w_gps]...
        = vcest(res_gps,fitmodel,vcmmodel,var0_gps_est,pb_gps,pl_gps,nbb_gps,w_gps);
    fprintf(fid,'current and cumulative GPS postpriori variance of unit weight: %5.2f, %5.2f\n',var0_gps,var0_gps_est);
    
    %  var0_gps = sparse(var0_gps);
    %  vcmmodel = sparse(vcmmodel);
    
    delvce=0;
    ksqr_insar=0;
    if ninsar>0
        res_insar=res(1:nobsinsar);
    end
    
    if vce==0 || ninsar==0
        vcmmodel = var0_gps*vcmmodel; %reweight using estimated variance of unit weight
    else
        %unit variance for InSAR
        if ninsar>0
            [ksqr_insar,var0_insar,var0_insar_est,pb_insar,pl_insar,nbb_insar,w_insar]...
                = vcest(res_insar,fitmodel,vcmmodel,var0_insar_est,pb_insar,pl_insar,nbb_insar,w_insar);
            delvce = delvce+abs(var0_gps/var0_insar-1);
            fprintf(fid,'current and cumulative InSAR postpriori variance of unit weight: %5.2f, %5.2f\n',var0_insar,var0_insar_est);
        end
        fprintf(fid,'variance component estimation - interation=%3d, delvce=%4.2f\n',iter,delvce);
        
        iter=iter+1;
    end
end

%--------------------------------------------
% output RMS vs Roughness
%--------------------------------------------
%RMS GPS
if ngf>0
    rms_gps=sqrt(mean(res_gps.^2));
    i1=cumsum([gps.nsite].*[gps.ndim]);
    i0=[1,i1(1:ngf-1)+1];
    for i=1:ngf
        ires=reshape(res_gps(i0(i):i1(i)),[],gps(i).ndim);
        rms_gps_enu=sqrt(mean(ires.^2));
        if gps(i).ndim==2
            rms_gps_enu(3)=nan;
        end
        fprintf(fid,'RMS_GPS_F%d (E, N, U, Total): %f %f %f %f\n',i,rms_gps_enu,rms_gps);
    end
end

%RMS InSAR
if ninsar>0
    rms_insar=sqrt(mean(res_insar.^2));
    fprintf(fid,'RMS_InSAR: %f\n',rms_insar);
end

%RMS
rms_all=sqrt(mean(res(1:m).^2));
fprintf(fid,'RMS ALL: %f\n',rms_all);

%wrss & roughness
if nargout>2
    wrss=sqrt((ksqr_gps+ksqr_insar)/m);
    rough=sqrt(res(m+1:nobs)'*res(m+1:nobs)/nsm)/smf;
    fprintf(fid,'weighted rss: %f\n',wrss);
    fprintf(fid,'roughness: %f\n',rough);
end

fclose(fid);
