
load /example/ATF/output/smf-1.50/precrash.mat;
[insar]=loadlics1lk(insarpar);
save('insar_1lk','insar','-v7.3');
smfdir = '/example/ATF/output/smf-1.50/';
load(strcat(smfdir, 'fitmodel.mat'));
[insarfit2]=insarfwd2(insar,trim,fitmodel,invenu,smfdir,gps);
cd (smfdir)
save('insarfit2','insarfit2','-v7.3');
cd ../
