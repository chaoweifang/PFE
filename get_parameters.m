function [param_ini,param_s,param_ucm]=get_parameters(iniflag,p, wtype, version)
%iniflag: character {g,h,s} corresponds to which type of initialization
%strategy g:GMM h:hierachical clustering s:spectral clustering

%all parameter for three types of initialization
param_ini=struct('flag',iniflag, 'radius1',3, 'spgf1',5,'sigf1',2,...
    'radius2',5,'spgf2',2,'sigf2',5,...                        
    'reg',200,'preg',0.8,'lb',50,'scale',1,...
    'chs',4,'awsc',5);

%all parameters for segmentation
param_s=struct('p',p,'wtype',wtype, 'miu',4000, 'r',600, 'rex',10, ...
    'mult_Pb',0.7, 'sat_sPb',0.8, 'dthresh',5,'ic_gamma',0.12,... 
    'or1',3,'or2',3, 'rscale',50,'rlb',0.01, 'version',version);
if strcmp(param_s.wtype,'o')
    param_s.rscale=1;
    param_s.rlb=1e-6;
    param_s.sat_sPb = param_s.sat_sPb/500;
end

param_ucm=struct('thr',0.35, 'fq',11);   
