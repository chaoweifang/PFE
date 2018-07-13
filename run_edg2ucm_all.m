addpath('./external/bsr/')
addpath('./external/structured_forest/')
addpath('./datasets/')
addpath('./lib')

clear clc

database='VOCContext';
gt_set='test';
res_root = '~/Results/';
res_dir=fullfile(res_root,database,'DB-UCM');
paras =struct('thr',0.5,'fq',8);

edg2ucm_all(database, gt_set,res_dir,paras)