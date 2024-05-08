function [proc] = getinsarproc(procfile)
%======================================================
%function [proc] = getinsarproc(procfile)
%
% Get processing parameters for insar dataset
%                                                      
% INPUT:                                               
%   procfile: filename of processing parameters (optional)       
% OUTPUT:                                              
%   proc: structure of processing parameters                            
%                                                      
% Hua Wang, 27/05/2011                         
%======================================================

%---------------------------------
%   DEFAULTS
%---------------------------------
%set parameters from global pars
if exist('pars.mat','file')
  load pars
  proc.orbdegree = orbdegree;
  proc.atmdegree = atmdegree;
  %modify here for compatible use in the interseismic inversion program
  if exist('invenu','var')
    proc.invenu = invenu;
  end
%set defaults
else
  proc.orbdegree = 1;
  proc.atmdegree = 1;
  proc.invenu = [1 1 0];
end

proc.errormap=0;
proc.incfile = 'incidence.unw';
proc.demfile = 'dem.dat';

%---------------------------------
%  CONFIGURE
%---------------------------------

%reading parameters from procfile
if nargin==1 & exist(procfile,'file')
  [parmat]=readparfile(procfile);

  %orbdegree: 0 - constant offset (reference phase)
  %           1 - linear orbital correction
  %           2 - quadrantic orbital correction (default)
  proc.orbdegree = getpar('orbdegree:',parmat,'n',proc.orbdegree);

  %atmdegree: 0 - without atm correction
  %           1 - linear topo-correlated atm correction
  proc.atmdegree = getpar('atmdegree:',parmat,'n',proc.atmdegree);

  %e/n/u contribution to the insar los rate
  if isfield(proc,'invenu')
    proc.invenu(1) = getpar('inv_e:',parmat,'n',proc.invenu(1));
    proc.invenu(2) = getpar('inv_n:',parmat,'n',proc.invenu(2));
    proc.invenu(3) = getpar('inv_u:',parmat,'n',proc.invenu(3));
  end

  %constant 1-sigma uncertainty of the rate map
  proc.errormap = getpar('errormap:',parmat,'n',proc.errormap);

  %incidence filename
  proc.incfile = getpar('incfile:',parmat,'s',proc.incfile);

  %dem filename
  proc.demfile = getpar('demfile:',parmat,'s',proc.demfile);
end
