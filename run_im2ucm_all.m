addpath('./pfe/');
addpath('./external/suitesparse/');
addpath('./external/bsr/')
addpath('./external/vlfeat/')
addpath('./scripts')
addpath('./external/structured_forest/')
addpath('./datasets/')
addpath('./lib')

clear;clc;

[param_ini,param_s,param_ucm]=get_parameters('s',0.8,'w', 'journal');
database='BSDS500';
gt_set='test';

res_root = '~/Results/';
res_dir=fullfile(res_root,database,['DB' param_ini.flag  param_s.wtype 'PFEp' num2str(param_s.p) '_' param_s.version]);
im2ucmpnorm_all(database, gt_set,res_dir,param_s,param_ucm,param_ini);
